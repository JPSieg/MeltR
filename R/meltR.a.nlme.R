#'Experimental code: not ready for use
#'
#'@param data_frame data_frame containing absorbance melting data
#'@param blank the blank sample
#'@param NucAcid A vector containing the Nucleic acid type and the sequences you are fitting.
#'@param Mmodel The molecular model you want to fit. Options: "Monomolecular.2State", "Monomolecular.3State", "Heteroduplex.2State", "Homoduplex.2State".
#'@param Tmodel The thermodynamic model you want to fit. Options: "VantHoff". Default = "VantHoff".
#'@param concT the temperature used to calculate the NucAcid concentration. Default = 90.
#'@param Save_results What results to save. Options: "all" to save PDF plots and ".csv" formated tables of parameters, "some" to save ".csv" formated tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@param file_path Path to the directory you want to save results in.
#'@return A list of data frames containing parameters from the fits and data for ploting the results with ggplot2.
#' @export
meltR.A.nlme = function(data_frame,
                   blank,
                   NucAcid,
                   concT = 90,
                   Mmodel,
                   Tmodel = "VantHoff",
                   Save_results = "none",
                   file_prefix = "Fit",
                   file_path = getwd()) {
  ####List of molecular models to fit####
  Mmodel_names <- c("Monomolecular.2State",
                    "Monomolecular.3State",
                    "Heteroduplex.2State",
                    "Homoduplex.2State")
  Mmodels <- list(function(K){ (K/(1+K)) },
                  function(K1, K2, Ct){ 1/(1 + K1 + (K1*K2)) },
                  function(K, Ct){ ((2/(K*Ct)) + 2 - sqrt(((((2/(K*Ct)) + 2)^2) - 4)))/2 },
                  function(K, Ct){ ((1/(2*K*Ct)) + 2 - sqrt(((((1/(2*K*Ct)) + 2)^2) - 4)))/2 })
  names(Mmodels) <- Mmodel_names
  ####List of thermodynamics models to fit to####
  Tmodel_names <- c("VantHoff")
  Tmodels <- list(function(H, Tm, Temperature){(H/0.0019872)*((1/(Tm + 273.15)) - (1/(Temperature + 273.15)))})
  names(Tmodels) <- Tmodel_names
  ####Assemble the models####
  G_VH = function(H, S, Temperature){exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H))}
  if (Mmodel == "Monomolecular.2State"){
    if (Tmodel == "VantHoff"){
      Model = function(H, Tm, mED, bED, mESS, bESS, Temperature){
        K <- exp(Tmodels$VantHoff(H = H, Tm = Tm, Temperature = Temperature))
        f <- Mmodels$Monomolecular.2State(K = K)
        ED <- mED*Temperature + bED
        ESS <- mESS*Temperature + bESS
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      GModel = function(H, S, mED, bED, mESS, bESS, Sample, Temperature){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        f <- Mmodels$Monomolecular.2State(K = K)
        ED <- mED[Sample]*Temperature + bED[Sample]
        ESS <- mESS[Sample]*Temperature + bESS[Sample]
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      calcS = function(H, Tm){ (H/(273.15 + Tm)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm)))) }
      calcG = function(H, Tm){H - (310*((H/(273.15 + Tm))))}
      calcG.SE = function(H, Tm, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*(H/(273.15 + Tm)))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2) }
    }
  }
  if (Mmodel == "Monomolecular.3State"){
    if (Tmodel == "VantHoff"){
      Model = function(H1, Tm1, H2, Tm2, mED, bED, EI, mESS, bESS, Temperature, Ct){
        K1 <- exp(Tmodels$VantHoff(H = H1, Tm = Tm1, Temperature = Temperature))
        K2 <- exp(Tmodels$VantHoff(H = H2, Tm = Tm2, Temperature = Temperature))
        f <- Mmodels$Monomolecular.2State(K1 = K1, K2 = K2, Ct = Ct)
        ED <- mED*Temperature + bED
        ESS <- mESS*Temperature + bESS
        model <- (f*K1*K2*ED) + (f*K1*EI) + (f)*ESS
        return(model)
      }
      GModel = function(H1, Tm1, H2, Tm2, mED, bED, EI, mESS, bESS, Sample, Temperature, Ct){
        K1 <- exp(Tmodels$VantHoff(H = H1, Tm = Tm1, Temperature = Temperature))
        K2 <- exp(Tmodels$VantHoff(H = H2, Tm = Tm2, Temperature = Temperature))
        f <- Mmodels$Monomolecular.2State(K1 = K1, K2 = K2, Ct = Ct)
        ED <- mED[Sample]*Temperature + bED[Sample]
        ESS <- mESS[Sample]*Temperature + bESS[Sample]
        model <- (f*K1*K2*ED) + (f*K1*EI[Sample]) + (f)*ESS
        return(model)
      }
      calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))}
    }
  }
  if (Mmodel == "Heteroduplex.2State"){
    if (Tmodel == "VantHoff"){
      Model = function(H, Tm, mED, bED, mESS, bESS, Temperature, Ct){
        K <- exp(Tmodels$VantHoff(H = H, Tm = Tm, Temperature = Temperature) + log(4/Ct))
        f <- Mmodels$Heteroduplex.2State(K = K, Ct = Ct)
        ED <- mED*Temperature + bED
        ESS <- mESS*Temperature + bESS
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      TmModel = function(H, S, lnCt){
        ((0.0019872/H)*lnCt) + ((S - 0.0019872*log(4))/H)
      }
      GModel = function(H, S, mED, bED, mESS, bESS, Sample, Temperature, Ct){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        f <- Mmodels$Heteroduplex.2State(K = K, Ct = Ct)
        ED <- mED[Sample]*Temperature + bED[Sample]
        ESS <- mESS[Sample]*Temperature + bESS[Sample]
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(4/Ct)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm)))) }
      calcG = function(H, Tm, Ct){ H - (310.15*((H/(273.15 + Tm)) + (0.0019872*log(4/Ct)))) }
      calcG.SE = function(H, Tm, Ct, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*((H/(273.15 + Tm)) + (0.0019872*log(4/Ct))))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2) }
    }
  }
  if (Mmodel == "Homoduplex.2State"){
    if (Tmodel == "VantHoff"){
      Model = function(H, Tm, mED, bED, mESS, bESS, Temperature, Ct){
        K <- exp(Tmodels$VantHoff(H = H, Tm = Tm, Temperature = Temperature) + log(1/Ct))
        f <- Mmodels$Homoduplex.2State(K = K, Ct = Ct)
        ED <- mED*Temperature + bED
        ESS <- mESS*Temperature + bESS
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      TmModel = function(H, S, lnCt){
        ((0.0019872/H)*lnCt) + (S/H)
      }
      GModel = function(H, S, mED, bED, mESS, bESS, Sample, Temperature, Ct){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        f <- Mmodels$Homoduplex.2State(K = K, Ct = Ct)
        ED <- mED[Sample]*Temperature + bED[Sample]
        ESS <- mESS[Sample]*Temperature + bESS[Sample]
        model <- (f*ED) + (1-f)*ESS
        return(model)
      }
      calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(1/Ct)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))}
      calcG = function(H, Tm, Ct){ H - (310.15*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct)))) }
      calcG.SE = function(H, Tm, Ct, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct))))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2)}
    }
  }
  ####Subtract out the blank####
  samples <- {}
  for (i in c(1:length(unique(data_frame$Sample)))){
    samples[[i]] <- subset(data_frame, Sample == unique(data_frame$Sample)[i])
    for (j in c(1:length(samples[[i]]$Sample))){
      samples[[i]]$Absorbance[j] <- samples[[i]]$Absorbance[j] - (samples[[blank]]$Absorbance[j]*samples[[blank]]$Pathlength[j])
    }
  }
  k <- samples[[1]]
  for (i in c(2:length(samples))){
    k <- rbind(k, samples[[i]])
  }
  no.background <- subset(k, Sample != blank)
  ####Calculate extinction coefficients####
  RNA <- list(15340, 7600, 12160, 10210, 13650, 10670, 12790, 12140, 10670, 7520, 9390, 8370, 12920, 9190, 11430, 10960, 12520, 8900, 10400, 10110)
  names(RNA) <- c("Ap", "Cp", "Gp", "Up", "ApA", "ApC", "ApG", "ApU", "CpA", "CpC", "CpG", "CpU", "GpA", "GpC", "GpG", "GpU", "UpA", "UpC", "UpG", "UpU")
  DNA <- list(15340, 7600, 12160, 8700, 13650, 10670, 12790, 11420, 10670, 7520, 9390, 7660, 12920, 9190, 11430, 10220, 11780, 8150, 9700, 8610)
  names(DNA) <- c("Ap", "Cp", "Gp", "Tp", "ApA", "ApC", "ApG", "ApT", "CpA", "CpC", "CpG", "CpT", "GpA", "GpC", "GpG", "GpT", "TpA", "TpC", "TpG", "TpT")
  if (NucAcid[1] == "RNA"){
    b <- c()
    b <- strsplit(NucAcid[-which(NucAcid == NucAcid[1])], split = "")
    c <- {}
    d <- {}
    e <- c()
    for (i in c(1:length(b))){
      c[[i]] <- c(1:(length(b[[i]])-1))
      d[[i]] <- c(1:(length(b[[i]])))
      for (j in c(1:(length(b[[i]])-1))){
        c[[i]][j] <- RNA[[which(names(RNA) == paste(b[[i]][j], "p", b[[i]][j+1], sep = ""))]]
      }
      for (j in c(1:(length(b[[i]])))){
        d[[i]][j] <- RNA[[which(names(RNA) == paste(b[[i]][j], "p", sep = ""))]]
      }
      e[i] <- 2*sum(c[[i]]) - sum(d[[i]])
    }
    extcoef <- list()
    extcoef[[1]] <- sum(e)
    for (i in c(1:length(b))){
      extcoef[[i+1]] <- e[i]
    }
    names(extcoef) <- c("Total", NucAcid[c(2:length(NucAcid))])
  }
  if (NucAcid[1] == "DNA"){
    b <- c()
    b <- strsplit(NucAcid[-which(NucAcid == NucAcid[1])], split = "")
    c <- {}
    d <- {}
    e <- c()
    for (i in c(1:length(b))){
      c[[i]] <- c(1:(length(b[[i]])-1))
      d[[i]] <- c(1:(length(b[[i]])))
      for (j in c(1:(length(b[[i]])-1))){
        c[[i]][j] <- RNA[[which(names(DNA) == paste(b[[i]][j], "p", b[[i]][j+1], sep = ""))]]
      }
      for (j in c(1:(length(b[[i]])))){
        d[[i]][j] <- RNA[[which(names(DNA) == paste(b[[i]][j], "p", sep = ""))]]
      }
      e[i] <- 2*sum(c[[i]]) - sum(d[[i]])
    }
    extcoef <- list()
    extcoef[[1]] <- sum(e)
    for (i in c(1:length(b))){
      extcoef[[i+1]] <- e[i]
    }
    names(extcoef) <- c("Total", NucAcid[c(2:length(NucAcid))])
  }
  ####Calculate Ct for each curve####
  samples <- {}
  if (length(extcoef) == 3){
    ct <- c()
    for (i in c(1:length(unique(no.background$Sample)))){
      samples[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
      for (j in c(length(samples[[i]]$Sample):1)){
        if (samples[[i]]$Temperature[j] > concT){
          ct[i] <- (samples[[i]]$Absorbance[j]/(extcoef$Total*samples[[i]]$Pathlength[j]))
        }
        samples[[i]]$Ct <- ct[i]
      }
    }
  }
  if (length(extcoef) == 2){
    ct <- c()
    for (i in c(1:length(unique(no.background$Sample)))){
      samples[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
      for (j in c(length(samples[[i]]$Sample):1)){
        if (samples[[i]]$Temperature[j] > concT){
          ct[i] <- (samples[[i]]$Absorbance[j]/(extcoef$Total*samples[[i]]$Pathlength[j]))
        }
        samples[[i]]$Ct <- ct[i]
      }
    }
  }
  k <- samples[[1]]
  for (i in c(2:length(samples))){
    k <- rbind(k, samples[[i]])
  }
  no.background <- subset(k, Sample != blank)
  ####Calculate starting thermo parameters for nls####
  first.derive <- {}
  T0.5 <- c()
  T0.75 <- c()
  startH <- c()
  for (i in c(1:length(unique(no.background$Sample)))){
    first.derive[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
    y <- c()
    x <- c()
    for (j in c(10:length(first.derive[[i]]$Temperature))){
      y[j] <- mean(first.derive[[i]]$Absorbance[j:(j-9)])
      x[j] <- mean(first.derive[[i]]$Temperature[j:(j-9)])
    }
    a <- c()
    b <- c()
    for (j in c(8:length(first.derive[[i]]$Temperature))){
      a[j] <- (y[j] - y[j-7])/(x[j] - x[j-7])
      b[j] <- (x[j] + x[j-7])/2
    }
    T0.5[i] <- b[which.max(a)]
    T0.75[i] <- min(b[which(a <= 0.5*max(a, na.rm = TRUE))][which(b[which(a <= 0.5*max(a, na.rm = TRUE))] > b[which.max(a)])], na.rm = TRUE)
    if (Mmodel == "Monomolecular.2State"){
      startH[i] <- -0.0032/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }
    if (Mmodel == "Heteroduplex.2State"){
      startH[i] <- -0.007/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }
    if (Mmodel == "Homoduplex.2State"){
      startH[i] <- -0.0044/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }
  }
  ####Calculate starting baseline values for nls####
  a <-{}
  uppbl_fit <- {}
  lowbl_fit <- {}
  bl.start.plot <- {}
  for (i in c(1:length(unique(no.background$Sample)))){
    a[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
    upperbl <- data.frame(
      "Absorbance" = a[[i]]$Absorbance[which(a[[i]]$Absorbance >= quantile(a[[i]]$Absorbance, probs = 0.75))],
      "Temperature" = a[[i]]$Temperature[which(a[[i]]$Absorbance >= quantile(a[[i]]$Absorbance, probs = 0.75))]
    )
    lowbl <- data.frame(
      "Absorbance" = a[[i]]$Absorbance[which(a[[i]]$Absorbance <= quantile(a[[i]]$Absorbance, probs = 0.25))],
      "Temperature" = a[[i]]$Temperature[which(a[[i]]$Absorbance <= quantile(a[[i]]$Absorbance, probs = 0.25))]
    )
    uppbl_fit[[i]] <- lm(Absorbance ~ Temperature, data = upperbl)
    lowbl_fit[[i]] <- lm(Absorbance ~ Temperature, data = lowbl)
  }
  ####Method 1 fit each curve individually####
  a <-{}
  fit <- {}
  indvfits.H <- c()
  indvfits.S <- c()
  indvfits.G <- c()
  indvfits.Tm <- c()
  bED <- c()
  mED <- c()
  bSS <- c()
  mSS <- c()
  for (i in c(1:length(unique(no.background$Sample)))){
    tryCatch({
      a[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
      fitstart <- list(H = startH[i], Tm = T0.5[i],
                       mED = lowbl_fit[[i]]$coefficients[2], bED = lowbl_fit[[i]]$coefficients[1],
                       mESS = uppbl_fit[[i]]$coefficients[2], bESS = uppbl_fit[[i]]$coefficients[1])
      if (Mmodel == "Monomolecular.2State"){
        fit[[i]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature),
                        data = a[[i]],
                        start = fitstart,
                        trace = FALSE,
                        nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
      }
      if (Mmodel == "Heteroduplex.2State"){
        fit[[i]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature, Ct),
                        data = a[[i]],
                        start = fitstart,
                        trace = FALSE,
                        nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
      }
      if (Mmodel == "Homoduplex.2State"){
        fit[[i]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature, Ct),
                        data = a[[i]],
                        start = fitstart,
                        trace = FALSE,
                        nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
      }
      mED[i] <- coef(fit[[i]])[3]
      bED[i] <- coef(fit[[i]])[4]
      mSS[i] <- coef(fit[[i]])[5]
      bSS[i] <- coef(fit[[i]])[6]
      indvfits.H[i] <- coef(fit[[i]])[1]
      if (Mmodel == "Monomolecular.2State"){
        indvfits.S[i] <- calcS(coef(fit[[i]])[1], coef(fit[[i]])[2])
        indvfits.G[i] <- calcG(coef(fit[[i]])[1], coef(fit[[i]])[2])
      }
      if (Mmodel == "Heteroduplex.2State"){
        indvfits.S[i] <- calcS(coef(fit[[i]])[1], coef(fit[[i]])[2], a[[i]]$Ct[1])
        indvfits.G[i] <- calcG(coef(fit[[i]])[1], coef(fit[[i]])[2], a[[i]]$Ct[1])
      }
      if (Mmodel == "Homoduplex.2State"){
        indvfits.S[i] <- calcS(coef(fit[[i]])[1], coef(fit[[i]])[2], a[[i]]$Ct[1])
        indvfits.G[i] <- calcG(coef(fit[[i]])[1], coef(fit[[i]])[2], a[[i]]$Ct[1])
      }
      indvfits.Tm[i] <- coef(fit[[i]])[2]
    }, error = function(e){print(a[[i]]$Sample[1])
      #mED[i] <- NA
      #bED[i] <- NA
      #mSS[i] <- NA
      #bSS[i] <- NA
      #indvfits.H[i] <- NA
      #indvfits.S[i] <-NA
      #indvfits.G[i] <- NA
  })}
  indvfits <- data.frame("Sample" = unique(no.background$Sample),
                         "Ct" = unique(no.background$Ct),
                         "H" = round(indvfits.H, 1),
                         "S" = round(indvfits.S, 4),
                         "G" = round(indvfits.G, 1),
                         "Tm" = round(indvfits.Tm, 1))
  indvfits.mean <- list(round(mean(indvfits$H), 1), round(sd(indvfits$H), 1), round(1000*mean(indvfits$S), 1), round(1000*sd(indvfits$S), 1), round(mean(indvfits.G), 1), round(sd(indvfits.G), 1))
  names(indvfits.mean) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
  indvfits.mean <- data.frame(indvfits.mean)
  no.background$Ext <- no.background$Absorbance/(no.background$Pathlength*no.background$Ct)
  if (Save_results == "all"){
    pdf(paste(file_path, "/", file_prefix, "_method_1_raw_fit_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(no.background$Temperature, no.background$Absorbance,
         xlab = Temperature ~ (degree ~ C), ylab = "Absorbance",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    for (i in c(1:length(a))){
      lines(a[[i]]$Temperature, predict(fit[[i]]), col = "red")
    }
    dev.off()
    pdf(paste(file_path, "/", file_prefix, "_method_1_normalized_fit_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(no.background$Temperature, no.background$Ext,
         xlab = Temperature ~ (degree ~ C), ylab = "Absorbtivity (1/M*cm)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    for (i in c(1:length(a))){
      lines(a[[i]]$Temperature, (predict(fit[[i]])/(a[[i]]$Ct[1]*a[[i]]$Pathlength[1])), col = "red")
    }
    dev.off()
  }
  ####Method nlme to test experimental effects####
  head(no.background)

  ?nlme::nlme
  ?nlme::nlmeControl

  nlme.fit <- nlme::nlme(Absorbance ~ ((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2 - sqrt(((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2)^2) - 4)))/2)*(mED*Temperature + bED)) + ((1-(((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2 - sqrt(((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2)^2) - 4)))/2))*(mESS*Temperature + bESS)),
                         fixed = list(S ~ 1, H ~ 1, mED ~ 1, bED ~ 1, mESS ~ 1, bESS ~ 1),
                         random = mED + bED + mESS + bESS ~ 1 | Sample,
                         start = c(S = mean(indvfits.S), H = mean(indvfits.H), mED = mean(mED), bED = mean(bED), mESS = mean(mSS), bESS = mean(bSS)),
                         data = no.background,
                         verbose = TRUE,
                         control = nlme::nlmeControl(returnObject = TRUE))
  nlme.fit <- nlme::nlme(Absorbance ~ ((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2 - sqrt(((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2)^2) - 4)))/2)*(mED*Temperature + bED)) + ((1-(((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2 - sqrt(((((2/((exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H)))*Ct)) + 2)^2) - 4)))/2))*(mESS*Temperature + bESS)),
                         fixed = list(S ~ 1, H ~ 1, mED ~ 1, bED ~ 1, mESS ~ 1, bESS ~ 1),
                         random = H + S + mED + bED + mESS + bESS ~ 1 | ROX,
                         start = c(S = mean(indvfits.S), H = mean(indvfits.H), mED = mean(mED), bED = mean(bED), mESS = mean(mSS), bESS = mean(bSS)),
                         data = no.background,
                         verbose = TRUE,
                         control = nlme::nlmeControl(returnObject = TRUE))

  summary(nlme.fit)

  ####Assemble Results####
  if (Mmodel == "Monomolecular.2State"){
    comparison <- rbind(indvfits.mean, Gfit_summary)
    comparison <- cbind(data.frame("Method" = c("1 individual fits", "3 Global fit")), comparison)
    row.names(comparison) <- c(1:2)
  }
  if (Mmodel == "Heteroduplex.2State"){
    comparison <- rbind(indvfits.mean, Tm_vs_lnCt, Gfit_summary)
    comparison <- cbind(data.frame("Method" = c("1 individual fits", "2 Tm versus ln[Ct]", "3 Global fit")), comparison)
    row.names(comparison) <- c(1:3)
  }
  if (Mmodel == "Homoduplex.2State"){
    comparison <- rbind(indvfits.mean, Tm_vs_lnCt, Gfit_summary)
    comparison <- cbind(data.frame("Method" = c("1 individual fits", "2 Tm versus ln[Ct]", "3 Global fit")), comparison)
    row.names(comparison) <- c(1:3)
  }
  if (Save_results != "none"){
    write.table(comparison, paste(file_path, "/", file_prefix, "_summary.csv", sep = ""), sep = ",", row.names = FALSE)
  }
  if (Save_results != "none"){
    write.table(indvfits, paste(file_path, "/", file_prefix, "_method_1_individual_fits.csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print("Individual curves")
  print(indvfits)
  print("Summary")
  print(comparison)
  output <- list("Summary" = comparison,
                 "Method.1.indvfits" = indvfits)
  if (Mmodel != "Monomolecular.2State"){
    output$Method.2.fit <- Tm_vs_lnCt_fit
  }
  output$Method.3.fit <- gfit
  output <- output
}
