#'Fit fluorescence binding isotherms to obtain thermodynamic parameters
#'
#'Description
#'
#'@param data_frame data_frame containing fluorescence binding data
#'@param Tmodel The thermodynamic model you want to fit the data to
#'@param K_error Acceptable uncetainty in equillibrium constants for fitting to thermodynamic models. K_error = K standard error/K. Default = 0.25.
#'@param Start_K A K value to start non-linear regression. Default = 0.1.
#'@param Optomize_B_conc Deals with a fundamental experimental uncertainty in determination of A =fluorophore and B = Quencher concentrations in the experiment. If TRUE, meltR.f will optimize the quencher labeled strand based on the shape of low temperature isotherms with K values < 0.1.
#'@param Low_reading Used by the concentration optimization algorithm. The isotherm, or reading, that you want to use to optomize the concentration. Default = "auto" will use the lowest temperature reading.
#'@param A_conc Used by the concentration optimization algorithm. The experimentally determined concentration in microMolar of the A = fluorophore labeled strand in the stock solution used to prepare the experimental solutions.
#'@param B_conc Used by the concentration optimization algorithm. The experimentally determined concentration in microMolar of the B = quencher labeled strand in the stock solution used to prepare the experimental solutions.
#'@param A_dilution_factor Used by the concentration optimization algorithm. The dilution factors used to prepare the final samples from the A = fluorophore labeled strand. Final concentration in nanoMolar = Dilution_factor*concentration stock in micromolar. Default = 50.
#'@param B_dilution_factor Used by the concentration optimization algorithm. The dilution factors used to prepare the final samples from the B = quencher labeled strand. Final concentration in nanoMolar = Dilution_factor*concentration stock in micromolar. Default = 250, 200, 150, 100, 62.5, 50, 37.5, 25, 12.70492, 2.540984, 0.254098, 0.
#'@param low_K Used by the concentration optimization algorithm. A low K value in nanomolar, that is used to find an optimum ration between A and B strands in the experiment. Default = 0.01.
#'@param Save_results What results to save. Options: "all" to save PDF plots and ".csv" formated tables of parameters, "some" to save ".csv" formated tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@param file_path Path to the directory you want to save results in.
#'@return
#' @export
meltR.F = function(data_frame,
                   Tmodel = "VantHoff",
                   K_error = 0.25,
                   Start_K = 0.1,
                   Optimize_B_conc = TRUE,
                   Low_reading = "auto",
                   A_conc = 4,
                   B_conc = 4,
                   A_dilution_factor = c(50),
                   B_dilution_factor = c(250, 200, 150, 100, 62.5, 50, 37.5, 25, 12.70492, 2.540984, 0.254098, 0),
                   low_K = 0.01,
                   Save_results = "none",
                   file_prefix = "Fit",
                   file_path = getwd()){
  ####List of molecular models to fit####
  Mmodel <- function(K, A, B){ ((K+A+B)-(((K+A+B)^2)-(4*A*B))^(1/2))/(2*A) }
  ####List of thermodynamics models to fit to####
  Tmodel_names <- c("VantHoff")
  Tmodels <- list(function(H, S, Temperature){((1/((Temperature + 273.15)*0.0019872))*H) - (S/0.0019872)})
  names(Tmodels) <- Tmodel_names
  ####Assemble the models####
  Optimize = function(R, K, A, B, Fmax, Fmin){
    f <- Mmodel(K = K, A = A, B = B*R)
    model <- Fmax + (Fmin-Fmax)*f
    return(model)
  }
  if (Tmodel == "VantHoff"){
    Individual = Global = function(K, Fmax, Fmin, A, B){
      f <- Mmodel(K = K, A = A, B = B)
      model <- Fmax + (Fmin-Fmax)*f
      return(model)
    }
    Global = function(H, S, Fmax, Fmin, Reading, A, B, Temperature){
      K <- (10^9)*exp(Tmodels$VantHoff(H = H, S = S, Temperature = Temperature))
      f <- Mmodel(K = K, A = A, B = B)
      model <- Fmax[Reading] + (Fmin[Reading] - Fmax[Reading])*f
      return(model)
    }
    calcG = function(H, S){H - (310.15*S)}
    calcG.SE = function(SE.H, SE.S, covar){ sqrt((SE.H)^2 + (310.15*SE.S)^2 - 2*310.15*SE.H*SE.S*covar) }
  }
  ####Find starting Fmax and Fmin and optimize B conc####
  a <- {}
  Fmax.start <- c()
  Fmin.start <- c()
  for (i in c(1:length(unique(data_frame$Reading)))){
    a[[i]] <- subset(data_frame, Reading == unique(data_frame$Reading)[i])
    Fmax.start[i] <- mean(a[[i]]$Emission[ which(a[[i]]$Emission >= quantile(a[[i]]$Emission, probs =0.80)) ])
    Fmin.start[i] <- mean(a[[i]]$Emission[ which(a[[i]]$Emission <= quantile(a[[i]]$Emission, probs =0.2)) ])
  }
  ####Optimize the mole ratio of fluorophore labeled RNA to quencher labeled RNA####
  if (Optimize_B_conc == TRUE){
    if (Low_reading == "auto"){
      Low_reading <- data_frame$Reading[which.min(data_frame$Temperature)]
    }
    if (Low_reading != "auto"){
      Low_reading <- data_frame$Reading[Low_reading]
    }
    optomize_start <- list(R = 1, Fmax = Fmax.start[Low_reading], Fmin = Fmin.start[Low_reading])
    hockey_stick <- a[[Low_reading]]
    hockey_stick$low_K <- low_K
    hockey_stick_fit <- nls(Emission ~ Optimize(R, K = low_K, A, B, Fmax, Fmin),
                            low = 0, algorithm = 'port',
                            start = optomize_start, data = hockey_stick)
    R <- coef(hockey_stick_fit)[1]
  }
  ####Change the Q concentrations to an optimal F/Q ratio####
  if (Optimize_B_conc == TRUE){
    for (i in c(1:length(data_frame$Well))){
      for (j in c(1:length(B_dilution_factor))){
        if (round(data_frame$B[i]) == round(B_conc*B_dilution_factor[j])){ data_frame$B[i] <- B_conc*B_dilution_factor[j]*R}
      }
    }
  }
  ####Method 1 Fit individual isotherms####
  a <- {}
  fit <- {}
  temp <- c()
  K <- c()
  Fmax <- c()
  Fmin <- c()
  SE.K <- c()
  for (i in c(1:length(unique(data_frame$Reading)))){
    a[[i]] <- subset(data_frame, Reading == unique(data_frame$Reading)[i])
    Start_fit <- list(K = Start_K, Fmax = Fmax.start[i], Fmin = Fmin.start[i])
    tryCatch({
      fit[[i]] <- nls(Emission ~ Individual(K, Fmax, Fmin, A, B), low =0,
                      algorithm = 'port', start = Start_fit, data = a[[i]])
      temp[i] <- a[[i]]$Temperature[1]
      K[i] <- coef(fit[[i]])[1]
      SE.K[i] <- summary(fit[[i]])$coefficients[1,2]
      Fmax[i] <- coef(fit[[i]])[2]
      Fmin[i] <- coef(fit[[i]])[3]
    }, error = function(e){
      temp[i] <- a[[i]]$Temperature[1]
      K[i] <- NA
      SE.K[i] <- NA
      Fmax[i] <- NA
      Fmin[i] <- NA
    })}
  indvfits <- data.frame("Temperature" =  temp,
                         "K" = K,
                         "SE.K" = SE.K,
                         "Fmax" = Fmax,
                         "Fmin" = Fmin)
  indvfits$invT <- 1/(273.15 + indvfits$Temperature)
  indvfits$lnK <- log((10^-9)*indvfits$K)
  indvfits$SE.lnK <- indvfits$SE.K/indvfits$K
  vh_start = list(H = -70, S = -0.180)
  vh_plot_fit <- nls(lnK~Tmodels$VantHoff(H, S, Temperature),
                     start = vh_start,
                     data = indvfits[which(indvfits$SE.lnK <= K_error),])

  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_1_VH_plot.pdf", sep = ""),
        width = 4, height = 4, pointsize = 3)
    plot(indvfits$invT, indvfits$lnK,
         xlab = "1/Temperature" ~ (degree ~ C), ylab = "ln[ K ]")
    for(i in c(1:length(indvfits$Temperature))){
      arrows(y0 = indvfits$lnK[i] - indvfits$SE.lnK[i], y1 = indvfits$lnK[i] + indvfits$SE.lnK[i],
             x0 = indvfits$invT[i], x1 = indvfits$invT[i],
             length=0.02, angle=90, code = 3)
    }
    lines(indvfits$invT, Tmodels$VantHoff(H = coef(vh_plot_fit)[1], S = coef(vh_plot_fit)[2], Temperature = indvfits$Temperature), col = "red")
    dev.off()
  }
  VH_plot_summary <- list(round(coef(vh_plot_fit)[1], 1), round(summary(vh_plot_fit)$coefficients[1,2], 1),
                          round(1000*coef(vh_plot_fit)[2], 1), round(1000*summary(vh_plot_fit)$coefficients[2,2], 1),
                          round(calcG(coef(vh_plot_fit)[1], coef(vh_plot_fit)[2]), 1), round(calcG.SE(summary(vh_plot_fit)$coefficients[1,2], summary(vh_plot_fit)$coefficients[2,2], summary(vh_plot_fit)$cov.unscaled[1,2]*(summary(vh_plot_fit)$sigma^2)), 1))
  names(VH_plot_summary) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
  VH_plot_summary <- data.frame(VH_plot_summary)
  ####Method 2 Global fit####
  gfit_data <- data.frame("Helix" = c(), "Well" = c(), "Reading" = c(),
                          "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c())
  for (i in which(indvfits$SE.lnK <= K_error)){
    gfit_data <- rbind(gfit_data, subset(data_frame, Reading == i))
  }
  b <- data.frame("Helix" = c(), "Well" = c(), "Reading" = c(),
                  "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c())
  for (i in c(1:length(unique(gfit_data$Reading)))){
    a <- subset(gfit_data, Reading == unique(gfit_data$Reading)[i])
    a$Reading <- i
    b <- rbind(b, a)
  }
  gfit_data <- b
  gfit_start = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000, Fmax = Fmax[which(indvfits$SE.lnK <= K_error)], Fmin = Fmin[which(indvfits$SE.lnK <= K_error)])
  gfit <- nls(Emission ~ Global(H, S, Fmax, Fmin, Reading, A, B, Temperature),
              start = gfit_start,
              data = gfit_data)
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_2_Gfit_plot.pdf", sep = ""),
        width = 4, height = 4, pointsize = 3)
    plot(gfit_data$B, gfit_data$Emission,
         xlab = "[Quencher] (nM)", ylab = "Emission")
    for (i in c(1:length(unique(gfit_data$Reading)))){
      a <- subset(gfit_data, Reading == unique(gfit_data$Reading)[i])
      lines(c(1:ceiling(max(a$B))), Global(H = coef(gfit)[1],
                                           S = coef(gfit)[2],
                                           Fmax = coef(gfit)[i + 2],
                                           Fmin = coef(gfit)[i + 24],
                                           A = a$A[1],
                                           B = c(1:ceiling(max(a$B))),
                                           Temperature = a$Temperature[1]),
            col = "red")
    }
    dev.off()
  }
  Gfit_summary <- list(round(coef(gfit)[1], 1), round(summary(gfit)$coefficients[1,2], 1),
                       round(1000*coef(gfit)[2], 1), round(1000*summary(gfit)$coefficients[2,2], 1),
                       round(calcG(coef(gfit)[1], coef(gfit)[2]), 1), round(calcG.SE(summary(gfit)$coefficients[1,2], summary(gfit)$coefficients[2,2], summary(gfit)$cov.unscaled[1,2]*(summary(gfit)$sigma^2)), 1))
  names(Gfit_summary) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
  Gfit_summary <- data.frame(Gfit_summary)
  ####Save results####
  output <- rbind(VH_plot_summary, Gfit_summary)
  row.names(output) <- c(1:2)
  output <- cbind(data.frame("Method" =c("1 VH plot", "2 Global fit")), output)
  print(output)
  if (Save_results != "none"){
    write.table(output, paste(file_path, "/", file_prefix, "_summary.csv", sep = ""), sep = ",", row.names = FALSE)
  }
  output <- output
}
