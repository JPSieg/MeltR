#'Fit fluorescence binding isotherms to obtain thermodynamic parameters
#'
#'Automates the trivial but time-consuming tasks associated with non-linear regression.
#'Uses two non-linear regression methods to calculate thermodynamic parameters. Method 1
#'fits each isotherm individually then calculates thermodynamic parameters with a Van't Hoff
#'plot. Method 2 calculates thermodynamic parameters with a global fit, where H and S are constant
#'between isotherms and the Fmax and Fmin are allowed to float. Also includes an algorithm that
#'optimizes the mole ratio of fluorophore labeled strands to quencher labeled strands.
#'
#'@param df data frame containing fluorescence binding data
#'@param Kd_error_quantile Quantile for uncertainty in equilibrium constants for fitting to thermodynamic models. K_error = K standard error/K. Default = 0.25 or the 25% most accurate binding constants in the Kd_range (See "Kd_range").
#'@param Kd_range A custom acceptable nM Kd range to fit for your experiment. Options FALSE for no custom Kd_range or c(start, end) to set a Kd_range. Example: "Kd_range = c(5, 100)"
#'@param Start_K A Kd value to start non-linear regression. Default = 0.1.
#'@param vh_start A list of initial guesses for the helix folding enthalpy and entropy. By default vh_start = list(H = -70, S = -0.180)
#'@param Optimize_conc Deals with a fundamental experimental uncertainty in determination of A =fluorophore and B = Quencher concentrations in the experiment. If TRUE, meltR.f will optimize the concentration for the fluorophore labeled strand based on the shape of low temperature isotherms.
#'@param Low_reading Used by the concentration optimization algorithm. The isotherm, or reading, that you want to use to optimize the concentration. Default = "auto" will use the lowest temperature reading.
#'@param low_K Used by the concentration optimization algorithm. A low Kd value in nanomolar, that is used to find an optimum ration between A and B strands in the experiment. Default = FALSE to allow the low_K to float in the concentration optimization algorithm. It is best to use values between 0.1 and 10.
#'@param B.conc.Tm Only use quencher (or B strands) higher than this threshold in the 1/Tm versus lnCt fitting method, method 3
#'@param Save_results What results to save. Options: "all" to save PDF plots and ".csv" formatted tables of parameters, "some" to save ".csv" formatted tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@param file_path Path to the directory you want to save results in.
#'@param silent Set to TRUE to run in silent mode (which does not print results in the console). Good for running in loops. Default is TRUE.
#'@return A meltR.F fit object containing a list of data objects containing raw data, data, transformation, fit objects, and statistics from the fits for plotting, exporting, and advanced analysis.
#' \itemize{
#'   \item 1. VantHoff - A data frame containing the duplex formation energies.
#'   \item 2. K - A data frame containing the results from fitting each isotherm individually.
#'   \item 3. VH_method_1_fit - A nls object containing the fit obtained from the fit obtained from the Van't Hoff plot.
#'   \item 4. VH_method_2_fit - A nls object containing the global fit.
#'   \item 5. Raw_data - The raw data passed back out of MeltR.F with no modifications.
#'   \item 6. First_derivative - The first derivative of each sample. Useful for qualitative comparison of data between conditions..
#'   \item 7. Tms - The approximate Tm of each sample obtained from the maximum of the first derivative. Useful for qualitative comparison of data between solution conditions.
#'   \item 8. R - The mole ratio of fluorophore and quencher labeled RNA, used in the concentration optimization algorithm. The mole ratio “R” was labeled “X” in the theory section to avoid confusion with the gas constant.
#'   \item 9. Fractional_error_between_methods - The amount thermodynamic parameters vary between methods.
#' }
#' @export
meltR.F = function(df,
                   Kd_error_quantile = 0.25,
                   Kd_range = c(10, 1000),
                   Start_K = 0.1,
                   vh_start = list(H = -70, S = -0.180),
                   Optimize_conc = TRUE,
                   Low_reading = "auto",
                   low_K = FALSE,
                   B.conc.Tm = 250,
                   Save_results = "none",
                   file_prefix = "Fit",
                   file_path = getwd(),
                   silent = FALSE){
  ####Make sure Temperature, A, B, Reading, and Emission are numeric####

  df$Reading = as.numeric(df$Reading)
  df$Temperature = as.numeric(df$Temperature)
  df$A = as.numeric(df$A)
  df$B = as.numeric(df$B)
  df$Emission = as.numeric(df$Emission)

  ####Generate a model to fit to####

  Mmodel = function(K, A, B){ ((K+A+B)-(((K+A+B)^2)-(4*A*B))^(1/2))/(2*A) }
  Tmodel = function(H, S, Temperature){((1/((Temperature + 273.15)*0.0019872))*H) - (S/0.0019872)}

  Optimize = function(R, K, A, B, Fmax, Fmin){
    f = Mmodel(K = K, A = A, B = B*R)
    model = Fmax + (Fmin-Fmax)*f
    return(model)
  }

  Individual = Global = function(K, Fmax, Fmin, A, B){
      f = Mmodel(K = K, A = A, B = B)
      model = Fmax + (Fmin-Fmax)*f
      return(model)
  }

  Global = function(H, S, Fmax, Fmin, Reading, A, B, Temperature){
      K = (10^9)*exp(Tmodel(H = H, S = S, Temperature = Temperature))
      f = Mmodel(K = K, A = A, B = B)
      model = Fmax[Reading] + (Fmin[Reading] - Fmax[Reading])*f
      return(model)
  }

  calcG = function(H, S){H - (310.15*S)}
  calcG.SE = function(SE.H, SE.S, covar){ sqrt((SE.H)^2 + (310.15*SE.S)^2 - (2*310.15*covar)) }

  Tm.v.A.B = function(H, S, A, B){
    if (A > B){
      invT = -((0.0019872/H)*log(A-(0.5*B))) + (S/H) -((9*0.0019872*log(10))/H)
    }
    if (A == B){
      invT = ((0.0019872/H)*log(B)) + (S/H) - ((0.0019872/H)*log(2)) - ((9*0.0019872*log(10))/H)
    }
    if (A < B){
      invT = ((0.0019872/H)*log(B-(0.5*A))) + (S/H) -((9*0.0019872*log(10))/H)
    }
    invT = invT
  }

  ####Find starting Fmax and Fmin and optimize B conc####

  Fmax.start = c()
  Fmin.start = c()
  for (i in c(1:length(unique(df$Reading)))){
    df.reading = subset(df, Reading == unique(df$Reading)[i])
    Fmax.start[i] = mean(df.reading$Emission[ which(df.reading$Emission >= quantile(df.reading$Emission, probs =0.80, na.rm = TRUE)) ])
    Fmin.start[i] = mean(df.reading$Emission[ which(df.reading$Emission <= quantile(df.reading$Emission, probs =0.2, na.rm = TRUE)) ])
  }

  ####Optimize the mole ratio of fluorophore labeled RNA to quencher labeled RNA####

  if (Optimize_conc == TRUE){
    if (Low_reading == "auto"){
      Low_reading = df$Reading[which.min(df$Temperature)]
    }else{Low_reading = df$Reading[Low_reading]}

    if (low_K == FALSE){
      optomize_start = list(R = 1, K = 0.01, Fmax = Fmax.start[Low_reading], Fmin = Fmin.start[Low_reading])
      hockey_stick = subset(df, Reading == Low_reading)
      hockey_stick$low_K = low_K
      hockey_stick_fit = nls(Emission ~ Optimize(R, K, A, B, Fmax, Fmin),
                             #trace = TRUE,
                             low = 0, algorithm = 'port',
                             start = optomize_start, data = hockey_stick,
                             control = nls.control(warnOnly = TRUE))
      R = coef(hockey_stick_fit)[1]
      K.opt = coef(hockey_stick_fit)[2]
    }else{
      optomize_start = list(R = 1, Fmax = Fmax.start[Low_reading], Fmin = Fmin.start[Low_reading])
      hockey_stick = subset(df, Reading == Low_reading)
      hockey_stick$low_K = low_K
      hockey_stick_fit = nls(Emission ~ Optimize(R, K = low_K, A, B, Fmax, Fmin),
                             low = 0, algorithm = 'port',
                             start = optomize_start, data = hockey_stick,
                             control = nls.control(warnOnly = TRUE))
      R = coef(hockey_stick_fit)[1]
      K.opt = low_K
    }

  }

  #plot(hockey_stick$B, hockey_stick$Emission)
  #lines(hockey_stick$B, predict(hockey_stick_fit),col="red",lty=2,lwd=3)

  ####Change the Q concentrations to an optimal F/Q ratio####

  if (Optimize_conc == TRUE){
    for (i in c(1:length(df$Well))){
        df$A[i] = df$A[i]/R
    }
  }else{
    R = NA
    K.opt = NA
  }


  ####Tm analysis####

  df.no.0 = subset(df, df$B >= 0.25*df$A)

  list.df.no.0 = {}

  Well = c()
  A = c()
  B = c()
  Tm = c()

  for (i in 1:length(unique(df.no.0$Well))){
    df.well = subset(df.no.0, Well == unique(df.no.0$Well)[i])
    fit.Em = lm(Emission ~ poly(Temperature, 20, raw=TRUE),
                data = df.well)

    string.dE.dT = gsub("NA", "0", gsub(", ", " + ", toString(paste(0:20, "*", coef(fit.Em), "*x^", -1:19, sep =""))))

    dE.dT = function(x){eval(parse(text = string.dE.dT))}

    #plot(df.well$Temperature, df.well$Emission)
    #lines(df.well$Temperature, predict(fit.Em), col = "red")
    #lines(seq(min(df.well$Temperature) + 2, max(df.well$Temperature) - 2, length.out = 1000),
    #      80*dE.dT(seq(min(df.well$Temperature) + 2, max(df.well$Temperature) - 2, length.out = 1000)), col = "red")
    #lines(seq(min(df.well$Temperature), max(df.well$Temperature), length.out = 1000),
    #      80*dE.dT(seq(min(df.well$Temperature), max(df.well$Temperature), length.out = 1000)), col = "blue")
    #abline(v = Tm[i])

    Well[i] = df.well$Well[1]
    A[i] = df.well$A[1]
    B[i] = df.well$B[1]
    df.well$dA.dT = dE.dT(df.well$Temperature)
    list.df.no.0[[i]] = df.well
    Tm[i] = seq(min(df.well$Temperature) + 2, max(df.well$Temperature) - 2, length.out = 1000)[which.max(dE.dT(seq(min(df.well$Temperature) + 2, max(df.well$Temperature) - 2, length.out = 1000)))]


  } #End for loop

  Tm_summary = data.frame(Well, A, B, Tm)
  Tm_summary$invT = 1/(273.15 + Tm_summary$Tm)

  #first.deriv.data
  df.no.0 = list.df.no.0[[1]]
  if (length(list.df.no.0) > 1){
    for (i in 2:length(df.no.0)){
      df.no.0 = rbind(df.no.0, list.df.no.0[[i]])
    }
  }

  ####Method 1 Fit individual isotherms####

  fit = {}
  Reading = c()
  temp = c()
  K = c()
  Fmax = c()
  Fmin = c()
  SE.K = c()

  for (i in c(1:length(unique(df$Reading)))){
    df.reading= subset(df, Reading == unique(df$Reading)[i])
    Start_fit = list(K = Start_K, Fmax = Fmax.start[i], Fmin = Fmin.start[i])
    tryCatch({
      fit[[i]] = nls(Emission ~ Individual(K, Fmax, Fmin, A, B), low =0,
                      algorithm = 'port', start = Start_fit, data = df.reading)
      temp[i] = df.reading$Temperature[1]
      K[i] = coef(fit[[i]])[1]
      SE.K[i] = summary(fit[[i]])$coefficients[1,2]
      Fmax[i] = coef(fit[[i]])[2]
      Fmin[i] = coef(fit[[i]])[3]
      Reading[i] = unique(df$Reading)[i]
    }, error = function(e){
      temp[i] = df.reading$Temperature[1]
      K[i] = NA
      SE.K[i] = NA
      Fmax[i] = NA
      Fmin[i] = NA
      Reading[i] = unique(df$Reading)[i]
    })}
  indvfits = data.frame("Temperature" =  temp,
                        "Reading" = Reading,
                         "K" = K,
                         "SE.K" = SE.K,
                         "Fmax" = Fmax,
                         "Fmin" = Fmin)
  indvfits$invT = 1/(273.15 + indvfits$Temperature)
  indvfits$lnK = log((10^-9)*indvfits$K)
  indvfits$SE.lnK = indvfits$SE.K/indvfits$K

  indvfits2 = data.frame("Temperature" =  temp,
                        "Reading" = Reading,
                        "K" = K,
                        "SE.K" = SE.K,
                        "Fmax" = Fmax,
                        "Fmin" = Fmin)
  indvfits2$invT = 1/(273.15 + indvfits$Temperature)
  indvfits2$lnK = log((10^-9)*indvfits$K)
  indvfits2$SE.lnK = indvfits$SE.K/indvfits$K

  if (length(which(is.na(indvfits$K))) > 0){
    indvfits = indvfits[-which(is.na(indvfits$K)),]
  }

  indvfits.to.fit = indvfits[-c(which(indvfits$K <= Kd_range[1]), which(indvfits$K >= Kd_range[2])),]

  vh_start = vh_start

  K_error = quantile(indvfits.to.fit$SE.lnK, Kd_error_quantile, na.rm = TRUE)

  indvfits.to.fit = indvfits.to.fit[-which(indvfits.to.fit$SE.lnK >= K_error),]

  vh_plot_fit = nls(lnK~Tmodel(H, S, Temperature),
                    start = vh_start,
                    control = nls.control(warnOnly = TRUE),
                    data = indvfits.to.fit)
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_1_VH_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(indvfits$invT, indvfits$lnK,
         xlab = "1/Temperature (K)", ylab = "ln[ KD (M)]",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    for(i in c(1:length(indvfits$Temperature))){
      arrows(y0 = indvfits$lnK[i] - indvfits$SE.lnK[i], y1 = indvfits$lnK[i] + indvfits$SE.lnK[i],
             x0 = indvfits$invT[i], x1 = indvfits$invT[i],
             length=0.02, angle=90, code = 3)
    }
    lines(indvfits$invT, Tmodel(H = coef(vh_plot_fit)[1], S = coef(vh_plot_fit)[2], Temperature = indvfits$Temperature), col = "red")
    abline(h = log((10^-9)*Kd_range[1]), col = "blue")
    abline(h = log((10^-9)*Kd_range[2]), col = "orange")
    dev.off()
  }

  VH_plot_summary = list(coef(vh_plot_fit)[1], summary(vh_plot_fit)$coefficients[1,2],
                          1000*coef(vh_plot_fit)[2], 1000*summary(vh_plot_fit)$coefficients[2,2],
                          calcG(coef(vh_plot_fit)[1], coef(vh_plot_fit)[2]), calcG.SE(summary(vh_plot_fit)$coefficients[1,2], summary(vh_plot_fit)$coefficients[2,2], summary(vh_plot_fit)$cov.unscaled[1,2]*(summary(vh_plot_fit)$sigma^2)))
  names(VH_plot_summary) = c("H", "SE.H", "S", "SE.S", "G", "SE.G")
  VH_plot_summary = data.frame(VH_plot_summary)

  indvfits = indvfits.to.fit

  ####Method 2 Global fit####

  df.gfit = data.frame("Helix" = c(), "Well" = c(), "Reading" = c(),
                          "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c(), "K" = c(), "index" = c())

  #Filter data by K error

  Fmax.start = c()
  Fmin.start = c()

  for (i in 1:nrow(indvfits)){
      df.reading = subset(df, Reading == indvfits$Reading[i])
      df.reading$Reading = i
      Fmax.start[i] = indvfits$Fmax[i]
      Fmin.start[i] = indvfits$Fmin[i]
      df.gfit = rbind(df.gfit, df.reading)
  }

  gfit_start = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000,
                    Fmax = Fmax.start,
                    Fmin = Fmin.start)

  gfit = nls(Emission ~ Global(H, S, Fmax, Fmin, Reading, A, B, Temperature),
             start = gfit_start,
             data = df.gfit,
             control = nls.control(warnOnly = TRUE))

  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_2_Gfit_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(df.gfit$B, df.gfit$Emission,
         xlab = "[Quencher] (nM)", ylab = "Emission",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    for (i in c(1:length(unique(df.gfit$Reading)))){
      a = subset(df.gfit, Reading == unique(df.gfit$Reading)[i])
      lines(c(1:ceiling(max(a$B))), Global(H = coef(gfit)[1],
                                           S = coef(gfit)[2],
                                           Fmax = coef(gfit)[i + 2],
                                           Fmin = coef(gfit)[i + 2 + length(unique(df.gfit$Reading))],
                                           A = a$A[1],
                                           B = c(1:ceiling(max(a$B))),
                                           Temperature = a$Temperature[1]),
            col = "red")
    }
    dev.off()
  }
  Gfit_summary = list(coef(gfit)[1], summary(gfit)$coefficients[1,2],
                       1000*coef(gfit)[2], 1000*summary(gfit)$coefficients[2,2],
                       calcG(coef(gfit)[1], coef(gfit)[2]), calcG.SE(summary(gfit)$coefficients[1,2], summary(gfit)$coefficients[2,2], summary(gfit)$cov.unscaled[1,2]*(summary(gfit)$sigma^2)))
  names(Gfit_summary) = c("H", "SE.H", "S", "SE.S", "G", "SE.G")
  Gfit_summary = data.frame(Gfit_summary)

  ####Method 3 Tm versus LnCt####

  df.Tm = Tm_summary

  #df.Tm$Ct = (10^-9)*df.Tm$B - (10^-9)*0.5*df.Tm$A

  #df.Tm = subset(df.Tm, df.Tm$Ct > 0)

  #df.Tm$lnCt = log(df.Tm$Ct)

  #fit = nls(invT ~ (S/H) + (0.00198720425864083/H)*lnCt,
  #    data = subset(df.Tm, B >= B.conc.Tm),
  #    start = list(H = -70, S = -0.22))

  #H = coef(fit)[1]
  #SE.H = coef(summary(fit))[1,2]
  #S = coef(fit)[2]
  #SE.S = coef(summary(fit))[2,2]
  #G = H - (273.13 + 37)*S
  #SE.G = calcG.SE(SE.H, SE.S, summary(fit)$cov.unscaled[1,2]*(summary(fit)$sigma^2))

  #df.Tm.result = data.frame(H, SE.H, 1000*S, 1000*SE.S, G, SE.G)

  #colnames(df.Tm.result) = c("H", "SE.H", "S", "SE.S", "G", "SE.G")

  #if (Save_results == "all"){
  #  pdf(paste(file_path, "/", file_prefix, "_method_3_Tm_vs_lnCt_plot.pdf", sep = ""),
  #      width = 3, height = 3, pointsize = 0.25)
  #  plot(df.Tm$lnCt, df.Tm$invT,
  #       xlab = "ln[Ct (M)]", ylab = "1/Tm (1/K)",
   #      cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
   # abline(a = S/H, b = 0.00198720425864083/H, col = "red")
  #  dev.off()
  #}

  ####Save results####
  output = {}
  output[[1]] = rbind(VH_plot_summary, Gfit_summary)
  row.names(output[[1]]) = c(1:2)
  output[[1]] = cbind(data.frame("Method" =c("1 VH plot", "2 Global fit")), output)
  output[[1]]$K_error = c(K_error, K_error)
  output[[1]]$R = R
  output[[1]]$Kd.opt = K.opt
  if (silent == FALSE){
    print("Van't Hoff")
    print(paste("accurate Ks = ", length(indvfits[which(indvfits$SE.lnK <= K_error[1]),]$SE.lnK), sep = ""))
    print(output[[1]])
  }

  if (Save_results != "none"){
    write.table(output, paste(file_path, "/", file_prefix, "_VH_summary.csv", sep = ""), sep = ",", row.names = FALSE)
  }
  output[[2]] = data.frame("Temperature" = indvfits2$Temperature,
                            "K" = 1/((10^-9)*indvfits2$K),
                            "SE.K" = ((10^9)*indvfits2$SE.K)/(indvfits2$K^2),
                            "Fmax" = indvfits2$Fmax,
                            "Fmin" = indvfits2$Fmin)
  output[[3]] = vh_plot_fit
  output[[4]] = gfit
  output[[5]] = df
  output[[6]] = df.no.0
  output[[7]] = Tm_summary
  if (Optimize_conc == TRUE){
    output[[8]] = R
  }
  if (Optimize_conc == FALSE){
    output[[8]] = NA
  }
  output[[9]] = data.frame("H" = abs((range(output[[1]]$H)[1] - range(output[[1]]$H)[2])/mean(output[[1]]$H)),
                            "S" = abs((range(output[[1]]$S)[1] - range(output[[1]]$S)[2])/mean(output[[1]]$S)),
                            "G" = abs((range(output[[1]]$G)[1] - range(output[[1]]$G)[2])/mean(output[[1]]$G)))
  if (silent == FALSE){
    print("Fractional error between methods")
    print(output[[9]])
    print("dH and dG are reporterd in kcal/mol and dS is in cal/mol/K. Tms are in deg Celcius")
  }

  names(output) = c("VantHoff", #1
                     "K", #2
                     "VH_method_1_fit", #3
                     "VH_method_2_fit", #4
                     "Raw_data", #5
                     "First_derivative", #6
                     "Tms", #7
                     "R", #8
                     "Fractional_error_between_methods") #9
  output  = output
}
