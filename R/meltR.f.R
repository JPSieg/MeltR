#'Fit fluorescence binding isotherms to obtain thermodynamic parameters
#'
#'Automates the trivial but time-consuming tasks associated with non-linear regression.
#'Uses two non-linear regression methods to calculate thermodynamic parameters. Method 1
#'fits each isotherm individually then calculates thermodynamic parameters with a Van't Hoff
#'plot. Method 2 calculates thermodynamic parameters with a global fit, where H and S are constant
#'between isotherms and the Fmax and Fmin are allowed to float. Also includes an algorithm that
#'optimizes the mole ratio of fluorophore labeled strands to quencher labeled strands.
#'
#'@param data_frame data_frame containing fluorescence binding data
#'@param Tmodel The thermodynamic model you want to fit the data to. Options are "VantHoff" and "Kirchoff".
#'@param K_error Acceptable uncetainty in equillibrium constants for fitting to thermodynamic models. K_error = K standard error/K. Default = c(0.25, 0.25).
#'@param K-range A custom acceptable nM Kd range to fit for your experiment. Options FALSE for no custom K_range or c(start, end) to set a K_range. Example: "K_range = c(5, 100)"
#'@param Start_K A K value to start non-linear regression. Default = 0.1.
#'@param Optomize_B_conc Deals with a fundamental experimental uncertainty in determination of A =fluorophore and B = Quencher concentrations in the experiment. If TRUE, meltR.f will optimize the quencher labeled strand based on the shape of low temperature isotherms with K values < 0.1.
#'@param Low_reading Used by the concentration optimization algorithm. The isotherm, or reading, that you want to use to optomize the concentration. Default = "auto" will use the lowest temperature reading.
#'@param low_K Used by the concentration optimization algorithm. A low K value in nanomolar, that is used to find an optimum ration between A and B strands in the experiment. Default = 0.01.
#'@param B.conc.Tm Only use quencher (or B strands) higher than this threshold in the 1/Tm versus lnCt fitting method, method 3
#'@param Save_results What results to save. Options: "all" to save PDF plots and ".csv" formated tables of parameters, "some" to save ".csv" formated tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@param file_path Path to the directory you want to save results in.
#'@param Tm_smooth Number of Temperature points your want to smooth as you determine the Tm for various points. Default = 4.
#'@return
#' @export
meltR.F = function(data_frame,
                   Tmodel = "VantHoff",
                   K_range = FALSE,
                   K_error = c(0.5, 0.5),
                   Start_K = 0.1,
                   Optimize_conc = TRUE,
                   Low_reading = "auto",
                   low_K = 0.1,
                   B.conc.Tm = 200,
                   Save_results = "none",
                   file_prefix = "Fit",
                   file_path = getwd(),
                   Tm_smooth = 8){
  ####Make sure Temperature, A, B, Reading, and Emission are numeric####
  data_frame$Reading <- as.numeric(data_frame$Reading)
  data_frame$Temperature <- as.numeric(data_frame$Temperature)
  data_frame$A <- as.numeric(data_frame$A)
  data_frame$B <- as.numeric(data_frame$B)
  data_frame$Emission <- as.numeric(data_frame$Emission)
  ####List of molecular models to fit####
  Mmodel <- function(K, A, B){ ((K+A+B)-(((K+A+B)^2)-(4*A*B))^(1/2))/(2*A) }
  ####List of thermodynamics models to fit to####
  Tmodel_names <- c("VantHoff", "Kirchoff")
  Tmodels <- list(
    function(H, S, Temperature){((1/((Temperature + 273.15)*0.0019872))*H) - (S/0.0019872)},
    function(H, S, C, Temperature){((1/((Temperature + 273.15)*0.0019872))*H) - (S/0.0019872) - ((C/0.0019872)*((310.15/(Temperature + 273.15)) - 1 + log((Temperature + 273.15)/310.15)))})
  names(Tmodels) <- Tmodel_names
  ####Assemble the models####
  Optimize = function(R, K, A, B, Fmax, Fmin){
    f <- Mmodel(K = K, A = A, B = B*R)
    model <- Fmax + (Fmin-Fmax)*f
    return(model)
  }
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
  calcG.SE = function(SE.H, SE.S, covar){ sqrt((SE.H)^2 + (310.15*SE.S)^2 - (2*310.15*covar)) }
  Tm.v.A.B = function(H, S, A, B){
    if (A > B){
      invT <- -((0.0019872/H)*log(A-(0.5*B))) + (S/H) -((9*0.0019872*log(10))/H)
    }
    if (A == B){
      invT <- ((0.0019872/H)*log(B)) + (S/H) - ((0.0019872/H)*log(2)) - ((9*0.0019872*log(10))/H)
    }
    if (A < B){
      invT <- ((0.0019872/H)*log(B-(0.5*A))) + (S/H) -((9*0.0019872*log(10))/H)
    }
    invT <- invT
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
  if (Optimize_conc == TRUE){
    if (Low_reading == "auto"){
      Low_reading <- data_frame$Reading[which.min(data_frame$Temperature)]
    }else{Low_reading <- data_frame$Reading[Low_reading]}
    optomize_start <- list(R = 1, Fmax = Fmax.start[Low_reading], Fmin = Fmin.start[Low_reading])
    hockey_stick <- a[[Low_reading]]
    hockey_stick$low_K <- low_K
    hockey_stick_fit <- nls(Emission ~ Optimize(R, K = low_K, A, B, Fmax, Fmin),
                            low = 0, algorithm = 'port',
                            start = optomize_start, data = hockey_stick)
    R <- coef(hockey_stick_fit)[1]
  }
  ####Change the Q concentrations to an optimal F/Q ratio####
  if (Optimize_conc == TRUE){
    for (i in c(1:length(data_frame$Well))){
        data_frame$A[i] <- data_frame$A[i]/R
    }
  }
  ####Tm analysis####
  a <- data_frame[which(data_frame$B >= 0.25*data_frame$A),]
  A <- c()
  B <- c()
  Tm <- c()
  b <- {}
  d <- {}
  n <- Tm_smooth

  for (i in 1:length(unique(a$Well))){
    b <- subset(a, Well == unique(a$Well)[i])
    fit.Em = lm(Emission ~ poly(Temperature, 10, raw=TRUE),
                data = b)
    b$dE.dT = coef(fit.Em)[2] + 2*coef(fit.Em)[3]*b$Temperature +
      3*coef(fit.Em)[4]*b$Temperature^2 + 4*coef(fit.Em)[5]*b$Temperature^3 +
      5*coef(fit.Em)[6]*b$Temperature^4 + 6*coef(fit.Em)[7]*b$Temperature^5 +
      7*coef(fit.Em)[8]*b$Temperature^6 + 8*coef(fit.Em)[9]*b$Temperature^7 +
      9*coef(fit.Em)[10]*b$Temperature^8 + 10*coef(fit.Em)[11]*b$Temperature^9
    Ts = seq(25, 75, length.out = 500)
    dE.dT = coef(fit.Em)[2] + 2*coef(fit.Em)[3]*Ts +
      3*coef(fit.Em)[4]*Ts^2 + 4*coef(fit.Em)[5]*Ts^3 +
      5*coef(fit.Em)[6]*Ts^4 + 6*coef(fit.Em)[7]*Ts^5 +
      7*coef(fit.Em)[8]*Ts^6 + 8*coef(fit.Em)[9]*Ts^7 +
      9*coef(fit.Em)[10]*Ts^8 + 10*coef(fit.Em)[11]*Ts^9
    Tm[i] = Ts[which.max(dE.dT)]
    Well <- unique(a$Well)[i]
    A[i] <- b$A[1]
    B[i] <- b$B[1]
    x <- c()
    y <- c()
    for (j in n:length(b$Temperature)){
      y[j] <- (b$Emission[j] - b$Emission[(j-(n-1))])/(b$Temperature[j] - b$Temperature[j-(n-1)])
      x[j] <- mean(b$Temperature[(j - n + 1):j])
    }
    d[[i]] <- data.frame("Well" = Well,
                         "Temperature" = x,
                         "A" = b$A[1],
                         "B" = b$B[1],
                         "First.derivative" = y)
  }

  e <- d[[1]]
  for (i in 2:length(d)){
    e <- rbind(e, d[[i]])
  }
  Tm_data <- e
  Tm_summary <- data.frame("Well" = unique(a$Well),
                           "A" = A,
                           "B" = B,
                           "Tm" = Tm)
  Tm_summary$invT <- 1/(273.15 + Tm_summary$Tm)
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_first_derivative.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(Tm_data$Temperature, Tm_data$First.derivative,
         xlab = "Temperature" ~ (degree ~ C), ylab = "dF/dT",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    dev.off()
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
  if(length(K_range) > 1){
    indvfits.to.fit <- indvfits[which(dplyr::between(indvfits$K, K_range[1], K_range[2])), ]
  }else{
    indvfits.to.fit <- indvfits
  }
  vh_start = list(H = -70, S = -0.180)
  if (is.na(K_error[1]) == FALSE){
    if (length(which(indvfits.to.fit$SE.lnK <= K_error[1])) <= 5){
      print(paste("Only ", length(which(indvfits.to.fit$SE.lnK <= K_error[1])), " isotherms have an acceptable K error. Try increasing the tolerance threshold.", sep = ""))
    }
    vh_plot_fit <- nls(lnK~Tmodels$VantHoff(H, S, Temperature),
                       start = vh_start,
                       data = indvfits.to.fit[which(indvfits.to.fit$SE.lnK <= K_error[1]),])
  }
  if (is.na(K_error[1]) == TRUE){
    vh_plot_fit <- nls(lnK~Tmodels$VantHoff(H, S, Temperature),
                       start = vh_start,
                       data = indvfits.to.fit)
  }
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_1_VH_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(indvfits$invT, indvfits$lnK,
         xlab = "1/Temperature" ~ (degree ~ C), ylab = "ln[ K ]",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
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
                          "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c(), "K" = c(), "index" = c())
  #Filter data by K error
  if (is.na(K_error[2]) == FALSE){
    for (i in which(indvfits$SE.lnK <= K_error[1])){
      a <- subset(data_frame, Reading == i)
      a$K <- indvfits$K[i]
      a$index <- i
      gfit_data <- rbind(gfit_data, a)
    }
  }
  if (is.na(K_error[2]) == TRUE){
    gfit_data <- data_frame
  }
  #Filter data by K range
  if (length(K_range) > 1){
    gfit_data <- gfit_data[which(dplyr::between(gfit_data$K, K_range[1], K_range[2])), ]
  }
  b <- data.frame("Helix" = c(), "Well" = c(), "Reading" = c(),
                  "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c(), "K" = c(), "index" = c())
  for (i in c(1:length(unique(gfit_data$Reading)))){
    a <- subset(gfit_data, Reading == unique(gfit_data$Reading)[i])
    a$Reading <- i
    b <- rbind(b, a)
  }
  gfit_data <- b
  if (is.na(K_error[2]) == FALSE){
    gfit_start = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000, Fmax = Fmax[unique(gfit_data$index)], Fmin = Fmin[unique(gfit_data$index)])
  }
  if (length(K_range) > 1){
    gfit_start = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000, Fmax = Fmax[unique(gfit_data$index)], Fmin = Fmin[unique(gfit_data$index)])
  }
  if (is.na(K_error[2]) == TRUE){
    gfit_start = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000, Fmax = Fmax, Fmin = Fmin)
  }
  gfit <- nls(Emission ~ Global(H, S, Fmax, Fmin, Reading, A, B, Temperature),
              start = gfit_start,
              data = gfit_data)
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_2_Gfit_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(gfit_data$B, gfit_data$Emission,
         xlab = "[Quencher] (nM)", ylab = "Emission",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    for (i in c(1:length(unique(gfit_data$Reading)))){
      a <- subset(gfit_data, Reading == unique(gfit_data$Reading)[i])
      lines(c(1:ceiling(max(a$B))), Global(H = coef(gfit)[1],
                                           S = coef(gfit)[2],
                                           Fmax = coef(gfit)[i + 2],
                                           Fmin = coef(gfit)[i + 2 + length(unique(gfit_data$Reading))],
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

  ####Method 3 Tm versus LnCt####

  df.Tm = Tm_summary

  df.Tm$Ct = (10^-9)*df.Tm$B - (10^-9)*0.5*df.Tm$A

  df.Tm$lnCt = log(df.Tm$Ct)

  fit = nls(invT ~ (S/H) + (0.00198720425864083/H)*lnCt,
      data = subset(df.Tm, B >= B.conc.Tm),
      start = list(H = -70, S = -0.22))

  H = coef(fit)[1]
  SE.H = coef(summary(fit))[1,2]
  S = coef(fit)[2]
  SE.S = coef(summary(fit))[2,2]
  G = H - (273.13 + 37)*S
  SE.G = calcG.SE(SE.H, SE.S, summary(fit)$cov.unscaled[1,2]*(summary(fit)$sigma^2))

  df.Tm.result = data.frame(H, SE.H, 1000*S, 1000*SE.S, G, SE.G)

  colnames(df.Tm.result) = c("H", "SE.H", "S", "SE.S", "G", "SE.G")

  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_3_Tm_vs_lnCt_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(df.Tm$lnCt, df.Tm$invT,
         xlab = "ln[Ct (M)]", ylab = "1/Tm (1/K)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    abline(a = S/H, b = 0.00198720425864083/H, col = "red")
    dev.off()
  }

  ####Fit to a Kirchoff thermodynamic model Method 1####
  if (Tmodel == "Kirchoff"){
    KC_start = list(H = Gfit_summary$H, S = Gfit_summary$S/1000, C = 0)
    KC_plot_fit <- nls(lnK ~ Tmodels$Kirchoff(H, S, C, Temperature),
                       start = KC_start,
                       data = indvfits[which(indvfits$SE.lnK <= K_error[1]),])
    summary(KC_plot_fit)
    if (Save_results == "all"){
      pdf(paste(file_path, "/", file_prefix, "_method_2_KC_plot.pdf", sep = ""),
          width = 3, height = 3, pointsize = 0.25)
      plot(indvfits$invT, indvfits$lnK,
           xlab = "1/Temperature" ~ (degree ~ C), ylab = "ln[ K ]",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
      for(i in c(1:length(indvfits$Temperature))){
        arrows(y0 = indvfits$lnK[i] - indvfits$SE.lnK[i], y1 = indvfits$lnK[i] + indvfits$SE.lnK[i],
               x0 = indvfits$invT[i], x1 = indvfits$invT[i],
               length=0.02, angle=90, code = 3)
      }
      lines(indvfits$invT, Tmodels$Kirchoff(H = coef(KC_plot_fit)[1], S = coef(KC_plot_fit)[2], C = coef(KC_plot_fit)[3], Temperature = indvfits$Temperature), col = "blue")
      lines(indvfits$invT, Tmodels$VantHoff(H = coef(vh_plot_fit)[1], S = coef(vh_plot_fit)[2], Temperature = indvfits$Temperature), col = "red")
      dev.off()
    }
    KC_plot_summary <- list(round(coef(KC_plot_fit)[1], 1), round(summary(KC_plot_fit)$coefficients[1,2], 1),
                            round(1000*coef(KC_plot_fit)[2], 1), round(1000*summary(KC_plot_fit)$coefficients[2,2], 1),
                            round(1000*coef(KC_plot_fit)[3], 1), round(1000*summary(KC_plot_fit)$coefficients[3,2], 1),
                            round(calcG(coef(KC_plot_fit)[1], coef(KC_plot_fit)[2]), 1), round(calcG.SE(summary(KC_plot_fit)$coefficients[1,2], summary(KC_plot_fit)$coefficients[2,2], summary(KC_plot_fit)$cov.unscaled[1,2]*(summary(KC_plot_fit)$sigma^2)), 1))
    names(KC_plot_summary) <- c("H", "SE.H", "S", "SE.S", "C", "SE.C", "G", "SE.G")
    KC_plot_summary <- data.frame(KC_plot_summary)
  }
  ####Fit to a Kirchoff thermodynamic model Method 2####
  if (Tmodel == "Kirchoff"){
    Global.KC = function(H, S, C, Fmax, Fmin, Reading, A, B, Temperature){
      K <- (10^9)*exp(Tmodels$Kirchoff(H = H, S = S, C = C, Temperature = Temperature))
      f <- Mmodel(K = K, A = A, B = B)
      model <- Fmax[Reading] + (Fmin[Reading] - Fmax[Reading])*f
      return(model)
    }
    gfit_data <- data.frame("Helix" = c(), "Well" = c(), "Reading" = c(),
                            "Temperature" = c(), "B" = c(), "A" = c(), "Emission" = c())
    for (i in which(indvfits$SE.lnK <= K_error[2])){
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
    gfit_start.KC = list(H = VH_plot_summary$H, S = VH_plot_summary$S/1000, C = 0, Fmax = Fmax[which(indvfits$SE.lnK <= K_error[2])], Fmin = Fmin[which(indvfits$SE.lnK <= K_error[2])])
    gfit.KC <- nls(Emission ~ Global.KC(H, S, C, Fmax, Fmin, Reading, A, B, Temperature),
                   start = gfit_start.KC,
                   data = gfit_data)
    summary(gfit.KC)
    if (Save_results == "all"){
      pdf(paste(file_path, "/", file_prefix, "_KC_method_2_Gfit_plot.pdf", sep = ""),
          width = 3, height = 3, pointsize = 0.25)
      plot(gfit_data$B, gfit_data$Emission,
           xlab = "[Quencher] (nM)", ylab = "Emission",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
      for (i in c(1:length(unique(gfit_data$Reading)))){
        a <- subset(gfit_data, Reading == unique(gfit_data$Reading)[i])
        lines(c(1:ceiling(max(a$B))), Global.KC(H = coef(gfit.KC)[1],
                                                S = coef(gfit.KC)[2],
                                                C = coef(gfit.KC)[3],
                                                Fmax = coef(gfit.KC)[i + 3],
                                                Fmin = coef(gfit.KC)[i + 3 + length(unique(gfit_data$Reading))],
                                                A = a$A[1],
                                                B = c(1:ceiling(max(a$B))),
                                                Temperature = a$Temperature[1]),
              col = "red")
      }
      dev.off()
    }
    Gfit.KC_summary <- list(round(coef(gfit.KC)[1], 1), round(summary(gfit.KC)$coefficients[1,2], 1),
                         round(1000*coef(gfit.KC)[2], 1), round(1000*summary(gfit.KC)$coefficients[2,2], 1),
                         round(1000*coef(gfit.KC)[3], 1), round(1000*summary(gfit)$coefficients[3,2], 1),
                         round(calcG(coef(gfit.KC)[1], coef(gfit.KC)[2]), 1), 1) #round(calcG.SE(summary(gfit.KC)$coefficients[1,2], summary(gfit.KC)$coefficients[2,2], summary(gfit.KC)$cov.unscaled[1,2]*(summary(gfit.KC)$sigma^2)), 1))
    names(Gfit.KC_summary) <- c("H", "SE.H", "S", "SE.S", "C", "SE.C", "G", "SE.G")
    Gfit.KC_summary <- data.frame(Gfit.KC_summary)
  }
  ####Method 3 1/Tm versus B analysis####
  ####Save results####
  output <- {}
  output[[1]] <- rbind(VH_plot_summary, Gfit_summary, df.Tm.result)
  row.names(output[[1]]) <- c(1:3)
  output[[1]] <- cbind(data.frame("Method" =c("1 VH plot", "2 Global fit", "3 1/Tm vs lnCT")), output)
  print("Van't Hoff")
  print(paste("accurate Ks = ", length(indvfits[which(indvfits$SE.lnK <= K_error[1]),]$SE.lnK), sep = ""))
  print(output[[1]])
  if (Save_results != "none"){
    write.table(output, paste(file_path, "/", file_prefix, "_VH_summary.csv", sep = ""), sep = ",", row.names = FALSE)
  }
  output[[2]] <- data.frame("Temperature" = indvfits$Temperature,
                            "K" = 1/((10^-9)*indvfits$K),
                            "SE.K" = ((10^9)*indvfits$SE.K)/(indvfits$K^2),
                            "Fmax" = indvfits$Fmax,
                            "Fmin" = indvfits$Fmin)
  output[[3]] <- vh_plot_fit
  output[[4]] <- gfit
  output[[5]] <- data_frame
  output[[6]] <- Tm_data
  output[[7]] <- Tm_summary
  if (Optimize_conc == TRUE){
    output[[8]] <- R
  }
  if (Optimize_conc == FALSE){
    output[[8]] <- NA
  }
  output[[9]] <- data.frame("H" = abs((range(output[[1]]$H)[1] - range(output[[1]]$H)[2])/mean(output[[1]]$H)),
                            "S" = abs((range(output[[1]]$S)[1] - range(output[[1]]$S)[2])/mean(output[[1]]$S)),
                            "G" = abs((range(output[[1]]$G)[1] - range(output[[1]]$G)[2])/mean(output[[1]]$G)))
  print("Fractional error between methods")
  print(output[[9]])
  names(output) <- c("VantHoff",
                     "K",
                     "VH_method_1_fit",
                     "VH_method_2_fit",
                     "Raw_data",
                     "First_derivative",
                     "Tms",
                     "R",
                     "Fractional_error_between_methods")
  if (Tmodel == "Kirchoff"){
    output[[10]] <- rbind(KC_plot_summary, Gfit.KC_summary)
    row.names(output[[10]]) <- c(1:2)
    output[[10]] <- cbind(data.frame("Method" =c("1 KC plot", "2 Global fit")), output[[5]])
    print("Kirchoff")
    print(paste("accurate Ks = ", length(indvfits[which(indvfits$SE.lnK <= K_error[2]),]$SE.lnK), sep = ""))
    print(output[[10]])
    output[[11]] <- KC_plot_fit
    output[[12]] <- gfit.KC
    names(output) <- c("VantHoff",
                       "K",
                       "VH_method_1_fit",
                       "VH_method_2_fit",
                       "Raw_data",
                       "First_derivative",
                       "Tms",
                       "R",
                       "Fractional_error_between_methods",
                       "Kirchoff",
                       "KC_method_1_fit",
                       "KC_method_2_fit")
  }
  output <- output
}
