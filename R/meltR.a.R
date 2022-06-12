#'Fit absorbance melting curves to obtain thermodynamic parameters
#'
#'Automates the trivial but time-consuming tasks associated with non-linear regression.
#'Calculates extinction coefficients, subtracts out the baseline buffer readings, and
#'calculates the strand concentration, Ct, in each sample. Then uses three non-linear regression
#'methods to calculate thermodynamic parameters. Method 1 fits each melting curve individually
#'then reports the average H and S from all of the curves. Method 2 calculates the Tm for each
#'melting curve, and calculates thermodynamic parameters by fitting the relationship between Tm
#'and Ct. Method 3 calculates thermodynamic parameters with a global fit, where H and S are constant
#'between isotherms and the baslines are allowed to float. Also includes an algorithm that
#'optimizes the mole ratio of fluorophore labeled strands to quencher labeled strands.
#'
#'@param data_frame data_frame containing absorbance melting data
#'@param blank The blank sample for background subtraction, or a list of blanks to apply to different samples for background subtraction. "none" to turn off background subtraction. If there is a single blank in the data set, The identity of the blank, for example, blank = 1 or blanke = "blank". If there are multiple blanks in the data, blank = list(c("blank 1", "Sample 1"), c("blank 2", "sample 2")) and so on. Sample identifiers should be what they are in the data frame. If you need to figure out what the sample identifiers are, run unique(df$Sample), where df is the name of the R data frame you are using, in your R console.
#'@param NucAcid A vector containing the Nucleic acid type and the sequences you are fitting.
#'@param Mmodel The molecular model you want to fit. Options: "Monomolecular.2State", "Monomolecular.3State", "Heteroduplex.2State", "Homoduplex.2State".
#'@param Tmodel The thermodynamic model you want to fit. Options: "VantHoff". Default = "VantHoff".
#'@param concT The temperature used to calculate the NucAcid concentration. Default = 90.
#'@param fitTs Option to only fit certain temperature ranges for melting curves. Either a vector or a list. If this is set to a vector, meltR.A will only fit temperatures in this range for all melting curves Example = c(17, 75). If set to a list of vectors, meltR.A will change what values are fit for each curve. Example, list(c(0,100), c(17,75), .... , c(0,100)). The length of this list has to be the equal to the number of samples that will be fit.
#'@param methods what methods do you want to use to fit data. Default = c(TRUE, TRUE, TRUE). Can be true or false. Note, method 1 must be set to TRUE or the subsequent steps will not work.
#'@param Save_results What results to save. Options: "all" to save PDF plots and ".csv" formated tables of parameters, "some" to save ".csv" formated tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@param file_path Path to the directory you want to save results in.
#'@return A list of data frames containing parameters from the fits and data for ploting the results with ggplot2.
#' @export
meltR.A = function(data_frame,
                   blank = "none",
                   NucAcid,
                   concT = 90,
                   fitTs = NULL,
                   methods = c(TRUE, TRUE, TRUE),
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
                  function(K, Ct){ ( (2/(K*Ct)) + 2 - sqrt(((((2/(K*Ct)) + 2)^2) - 4)))/2 },
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

  if (is.list(blank)){ #routine for multiple blanks in a data set

    df.samples <- {}
    df.blank <- {}

    for (i in c(1:length(blank))){
      df.samples[[i]] = subset(data_frame, Sample == blank[[i]][1])
      df.blank[[i]] = subset(data_frame, Sample == blank[[i]][2])
      df.samples[[i]]$Absorbance = df.samples[[i]]$Absorbance - df.blank[[i]]$Absorbance
    }
    k <- df.samples[[1]]
    if (length(blank) > 1){
      for (i in c(2:length(df.samples))){
        k <- rbind(k, df.samples[[i]])
      }
    }

    no.background = k

  }else{
    if (blank == "none"){ #routine for no blank in the data set

      no.background = data_frame

    }else{  #routine for a single blank in the data set

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

    }
  }


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

  ####Remove values not in fitTs####

  if (is.null(fitTs) == FALSE){
    if (is.list(fitTs) == FALSE){
      ranges <- rep(list(fitTs), length(unique(no.background$Sample)))
    }else{
      ranges <- fitTs
    }
    a <- {}
    for (i in c(1:length(unique(no.background$Sample)))){
      a[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])

      if(length(which(a[[i]]$Temperature > ranges[[i]][2])) > 0){
        a[[i]] <- a[[i]][ -which(a[[i]]$Temperature > ranges[[i]][2]) ,]
      }

      if(length(which(a[[i]]$Temperature < ranges[[i]][1])) > 0){
        a[[i]] <- a[[i]][ -which(a[[i]]$Temperature < ranges[[i]][1]) ,]
      }



    }
    b <- a[[1]]
    if (length(a) > 1){
      for (i in c(2:length(unique(no.background$Sample)))){
        b <- rbind(b, a[[i]])
      }
    }
    no.background <- b
  }

  ####Calculate Ct for each curve####
  samples <- {}
  if (length(extcoef) == 3){
    ct <- c()
    for (i in c(1:length(unique(no.background$Sample)))){
      samples[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
      df.raw = subset(data_frame, Sample == unique(no.background$Sample)[i])
      for (j in c(length(df.raw$Sample):1)){
        if (df.raw$Temperature[j] > concT){
          ct[i] <- (df.raw$Absorbance[j]/(extcoef$Total*df.raw$Pathlength[j]))
        }
        samples[[i]]$Ct <- ct[i]
      }
    }
  }
  if (length(extcoef) == 2){
    ct <- c()
    for (i in c(1:length(unique(no.background$Sample)))){
      samples[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
      df.raw = subset(data_frame, Sample == unique(no.background$Sample)[i])
      for (j in c(length(df.raw$Sample):1)){
        if (df.raw$Temperature[j] > concT){
          ct[i] <- (df.raw$Absorbance[j]/(extcoef$Total*df.raw$Pathlength[j]))
        }
        samples[[i]]$Ct <- ct[i]
      }
    }
  }
  k <- samples[[1]]
  if (length(samples) > 1){
    for (i in c(2:length(samples))){
      k <- rbind(k, samples[[i]])
    }
  }
  no.background <- subset(k, Sample != blank)
  ####Calculate starting thermo parameters for nls####
  first.derive <- {}
  T0.5 <- c()
  T0.75 <- c()
  startH <- c()
  for (i in c(1:length(unique(no.background$Sample)))){
    first.derive[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i])
    fit.Em = lm(Absorbance ~ poly(Temperature, 20, raw=TRUE),
                data = first.derive[[i]])

    string.dE.dT = gsub("NA", "0", gsub(", ", " + ", toString(paste(0:20, "*", coef(fit.Em), "*x^", -1:19, sep =""))))
    string.dE.dT2 = gsub("NA", "0", gsub(", ", " + ", toString(paste(0:20, "*", c(0,0:19),"*", coef(fit.Em), "*x^", -2:18, sep =""))))

    dE.dT = function(x){eval(parse(text = string.dE.dT))}
    dE.dT2 = function(x){eval(parse(text = string.dE.dT2))}

    T0.5[i] = seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000)[which.max(dE.dT(seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000)))]
    T0.75[i] = seq(min(first.derive[[i]]$Temperature) + 10, max(first.derive[[i]]$Temperature) - 10, length.out = 1000)[which.min(dE.dT2(seq(min(first.derive[[i]]$Temperature) + 10, max(first.derive[[i]]$Temperature) - 10, length.out = 1000)))]

    if (Mmodel == "Monomolecular.2State"){
      startH[i] <- -0.0032/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }
    if (Mmodel == "Heteroduplex.2State"){
      startH[i] <- -0.007/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }
    if (Mmodel == "Homoduplex.2State"){
      startH[i] <- -0.0044/((1/(273.15 + T0.5[i]))-(1/(273.15 + T0.75[i])))
    }

    #svg("~/Desktop/Melt_curve_derivative.svg")

    #plot(first.derive[[i]]$Temperature, first.derive[[i]]$Absorbance,
    #     xlab = Temperature ~ (degree ~ C), ylab = "Absorbance")
    #lines(first.derive[[i]]$Temperature, predict(fit.Em), col = "red")
    #lines(seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000),
    #  1.4 +10*dE.dT(seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000)), col = "red")
    #lines(seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000),
    #      1.475+10*dE.dT2(seq(min(first.derive[[i]]$Temperature) + 2, max(first.derive[[i]]$Temperature) - 2, length.out = 1000)), col = "blue")
    #abline(v = T0.75[i])
    #abline(v = T0.5[i])

    #dev.off()

    first.derive[[i]]$dA.dT = dE.dT(first.derive[[i]]$Temperature) + dE.dT(first.derive[[i]]$Temperature)*((first.derive[[i]]$Absorbance - predict(fit.Em))/first.derive[[i]]$Absorbance)
    first.derive[[i]]$dA.dT2 = dE.dT2(first.derive[[i]]$Temperature) + dE.dT2(first.derive[[i]]$Temperature)*((first.derive[[i]]$Absorbance - predict(fit.Em))/first.derive[[i]]$Absorbance)

    #plot(first.derive[[i]]$Temperature, first.derive[[i]]$dA.dT)
    #plot(first.derive[[i]]$Temperature, first.derive[[i]]$dA.dT2)

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
  if (methods[1] == TRUE){
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
      }, error = function(e){print(a[[i]]$Sample[1])})
    }
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
  }
  if (methods[1] == FALSE){
    indvfits <- data.frame("Sample" = NA,
                           "Ct" = NA,
                           "H" = NA,
                           "S" = NA,
                           "G" = NA,
                           "Tm" = NA)
    indvfits.mean <- list(NA, NA, NA, NA, NA, NA)
    names(indvfits.mean) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
    indvfits.mean <- data.frame(indvfits.mean)
  }
  ####Method 2 Tm vs Ct####
  if (methods[2] == TRUE){
    a <- {}
    Tm_range <- {}
    Tm_fit <- {}
    Tm <- c()
    lnCt <- c()
    for (i in c(1:length(unique(no.background$Sample)))){
      a[[i]] <- subset(no.background, Sample == unique(no.background$Sample)[i]) #Pull out data
      #calculate f at each temp
      a[[i]]$f <- (a[[i]]$Absorbance - (coef(fit[[i]])[5]*a[[i]]$Temperature + coef(fit[[i]])[6]))/((coef(fit[[i]])[3]*a[[i]]$Temperature + coef(fit[[i]])[4]) - (coef(fit[[i]])[5]*a[[i]]$Temperature + coef(fit[[i]])[6]))
      Tm_range[[i]] <- data.frame( #Pull out data where f >= 0.4 and f <= 0.6
        "Temperature" = a[[i]]$Temperature[which(a[[i]]$f >= 0.4)][which(a[[i]]$f[which(a[[i]]$f >= 0.4)] <= 0.6)],
        "f" = a[[i]]$f[which(a[[i]]$f >= 0.4)][which(a[[i]]$f[which(a[[i]]$f >= 0.4)] <= 0.6)]
      )
      Tm_fit[[i]] <- lm(f~Temperature, data = Tm_range[[i]]) #fit data where f >= 0.4 and f <= 0.6
      Tm[i] <- (0.5 - coef(Tm_fit[[i]])[1])/coef(Tm_fit[[i]])[2]
      lnCt[i] <- log(a[[i]]$Ct[1])
    }
    Tm_data <- data.frame("lnCt" = lnCt, "Tm" = Tm)
    Tm_data$invT <- 1/(Tm_data$Tm + 273.15)
    if (Mmodel != "Monomolecular.2State"){
      if (Mmodel != "Monomolecular.3State"){
        Tm_vs_lnCt_fit <- nls(invT ~ TmModel(H, S, lnCt),
                              data = Tm_data,
                              start = list(H = -70, S = -0.17))
        Tm_vs_lnCt <- list(round(coef(Tm_vs_lnCt_fit)[1], 1), round(summary(Tm_vs_lnCt_fit)$coefficients[1,2], 1),
                           round(1000*coef(Tm_vs_lnCt_fit)[2], 1), round(1000*summary(Tm_vs_lnCt_fit)$coefficients[2,2], 1),
                           round(coef(Tm_vs_lnCt_fit)[1] - (310.15*coef(Tm_vs_lnCt_fit)[2]), 1), round(sqrt((summary(Tm_vs_lnCt_fit)$coefficients[1,2])^2 + (310*summary(Tm_vs_lnCt_fit)$coefficients[2,2])^2 - 2*310*summary(Tm_vs_lnCt_fit)$cov.unscaled[1,2]*(summary(Tm_vs_lnCt_fit)$sigma^2)),1))
        names(Tm_vs_lnCt) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
        Tm_vs_lnCt <- data.frame(Tm_vs_lnCt)
        if (Save_results == "all"){
          pdf(paste(file_path, "/", file_prefix, "_method_2_plot.pdf", sep = ""),
              width = 3, height = 3, pointsize = 0.25)
          plot(x = Tm_data$lnCt, y = Tm_data$invT,
               xlab = "ln[ Ct (M) ]", ylab = "1/[ Temperature (K) ]",
               cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
          lines(x = Tm_data$lnCt, predict(Tm_vs_lnCt_fit), col = "red")
          dev.off()
        }
      }
    }
    if (Mmodel == "Monomolecular.2State"){
      if (Mmodel == "Monomolecular.3State"){
        if (Save_results == "all"){
          pdf(paste(file_path, "/", file_prefix, "_method_2_plot.pdf", sep = ""),
              width = 3, height = 3, pointsize = 0.25)
          plot(x = Tm_data$lnCt, y = Tm_data$invT,
               xlab = "ln[ Ct (M) ]", ylab = "1/[ Temperature (K) ]",
               cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
          lines(Tm_data$lnCt, predict(lm(Tm_data$invT ~ Tm_data$lnCt)), col = "red")
          dev.off()
        }
      }}
  }
  if (methods[2] == FALSE){
    Tm_vs_lnCt <- list(NA, NA,
                       NA, NA,
                       NA, NA)
    names(Tm_vs_lnCt) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
    Tm_vs_lnCt <- data.frame(Tm_vs_lnCt)
  }
  ####Method 3 Global fitting####
  if (methods[3] == TRUE){
    b <- data.frame("Sample" = c(), "Pathlength" = c(), "Temperature" = c(),
                    "Absorbance" = c(), "Ct" = c(), "Ext" = c())
    for (i in c(1:length(unique(no.background$Sample)))){
      a <- subset(no.background, Sample == unique(no.background$Sample)[i])
      a$Sample <- i
      b <- rbind(b, a)
    }
    gfit_data <- b
    gfit_start = list(H = mean(indvfits$H), S = mean(indvfits$S), mED = mED, bED = bED, mESS = mSS, bESS = bSS)
    if (Mmodel == "Monomolecular.2State"){
      gfit <- nls(Absorbance ~ GModel(H, S, mED, bED, mESS, bESS, Sample, Temperature),
                  start = gfit_start,
                  data = gfit_data,
                  nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
    }
    if (Mmodel == "Heteroduplex.2State"){
      gfit <- nls(Absorbance ~ GModel(H, S, mED, bED, mESS, bESS, Sample, Temperature, Ct),
                  start = gfit_start,
                  data = gfit_data,
                  nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
    }
    if (Mmodel == "Homoduplex.2State"){
      gfit <- nls(Absorbance ~ GModel(H, S, mED, bED, mESS, bESS, Sample, Temperature, Ct),
                  start = gfit_start,
                  data = gfit_data,
                  nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
    }
    if (Save_results == "all"){
      pdf(paste(file_path, "/", file_prefix, "_method_3_Gfit_raw_plot.pdf", sep = ""),
          width = 3, height = 3, pointsize = 0.25)
      plot(gfit_data$Temperature, gfit_data$Absorbance,
           xlab = Temperature ~ (degree ~ C), ylab = "Absorbance",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
      for (i in c(1:length(unique(gfit_data$Sample)))){
        a <- subset(gfit_data, Sample == unique(gfit_data$Sample)[i])
        if (Mmodel == "Monomolecular.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature)))),
                col = "red")
        }
        if (Mmodel == "Heteroduplex.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature))),
                                                         Ct = a$Ct[1]),
                col = "red")
        }
        if (Mmodel == "Homoduplex.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature))),
                                                         Ct = a$Ct[1]),
                col = "red")
        }
      }
      dev.off()
      pdf(paste(file_path, "/", file_prefix, "_method_3_Gfit_normalized_plot.pdf", sep = ""),
          width = 3, height = 3, pointsize = 0.25)
      plot(gfit_data$Temperature, gfit_data$Ext,
           xlab = Temperature ~ (degree ~ C), ylab = "Absorbtivity (1/M*cm)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
      for (i in c(1:length(unique(gfit_data$Sample)))){
        a <- subset(gfit_data, Sample == unique(gfit_data$Sample)[i])
        if (Mmodel == "Monomolecular.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature))))/(a$Pathlength[1]*a$Ct[1]),
                col = "red")
        }
        if (Mmodel == "Heteroduplex.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature))),
                                                         Ct = a$Ct[1])/(a$Pathlength[1]*a$Ct[1]),
                col = "red")
        }
        if (Mmodel == "Homoduplex.2State"){
          lines(c(5:ceiling(max(a$Temperature))), GModel(H = coef(gfit)[1],
                                                         S = coef(gfit)[2],
                                                         mED = coef(gfit)[i + 2],
                                                         bED = coef(gfit)[i + 2 + length(unique(gfit_data$Sample))],
                                                         mESS = coef(gfit)[i + 2 + 2*length(unique(gfit_data$Sample))],
                                                         bESS = coef(gfit)[i + 2 + 3*length(unique(gfit_data$Sample))],
                                                         Temperature = c(5:ceiling(max(a$Temperature))),
                                                         Ct = a$Ct[1])/(a$Pathlength[1]*a$Ct[1]),
                col = "red")
        }

      }
      dev.off()
    }
    Gfit_summary <- list(round(coef(gfit)[1], 1), round(summary(gfit)$coefficients[1,2], 1),
                         round(1000*coef(gfit)[2], 1), round(1000*summary(gfit)$coefficients[2,2], 1),
                         round(coef(gfit)[1] - (310.15)*coef(gfit)[2], 1), round(sqrt((summary(gfit)$coefficients[1,2])^2 + (310.15*summary(gfit)$coefficients[2,2])^2 - 2*310.15*summary(gfit)$cov.unscaled[1,2]*(summary(gfit)$sigma^2)), 1))
    names(Gfit_summary) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
    Gfit_summary <- data.frame(Gfit_summary)
  }
  if (methods[3] == FALSE){
    Gfit_summary <- list(NA, NA,
                         NA, NA,
                         NA, NA)
    names(Gfit_summary) <- c("H", "SE.H", "S", "SE.S", "G", "SE.G")
    Tm_vs_lnCt <- data.frame(Gfit_summary)
  }
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

  ####Assemble data for custom derivative plots####

  df.deriv = first.derive[[1]]

  if (length(first.derive) > 1){
    for (i in 2:length(first.derive)){
      #print(i)
      df.deriv = rbind(df.deriv, first.derive[[i]])
    }
  }

  ####Assemble data for custom individual fit plots####

  df.method.1 = subset(no.background, Sample == unique(no.background$Sample)[1])
  df.method.1$Model = predict(fit[[1]])

  if (length(unique(no.background$Sample)) >1){
    for (i in 2:length(unique(no.background$Sample))){
      df.method.1.i = subset(no.background, Sample == unique(no.background$Sample)[i])
      df.method.1.i$Model = predict(fit[[i]])
      df.method.1 = rbind(df.method.1, df.method.1.i)
    }
  }


  ####Assemble data for custom 1/Tm versus lnCt plots####

  if (Mmodel != "Monomolecular.2State"){
    if (methods[2] == TRUE){
      Tm_data$Model = predict(Tm_vs_lnCt_fit)
    }
  }

  ####Assemble global fit data for custom fit plots####

  if (methods[3] == TRUE){
    gfit_data$Model = predict(gfit)
  }

  ####Assemble final output####

  print("Individual curves")
  print(indvfits)
  print("Summary")
  print(comparison)
  range <- data.frame("H" = abs((range(comparison$H)[1]-range(comparison$H)[2])/mean(comparison$H)),
                      "S" = abs((range(comparison$S)[1]-range(comparison$S)[2])/mean(comparison$S)),
                      "G" = abs((range(comparison$G)[1]-range(comparison$G)[2])/mean(comparison$G)))
  print("fractional error between methods")
  print(range)
  output <- list("Summary" = comparison,
                 "Method.1.indvfits" = indvfits,
                 "Range" = range,
                 "Derivatives.data" = df.deriv,
                 "Method.1.data" = df.method.1,
                 "Method.1.fit" = fit)
  if (Mmodel != "Monomolecular.2State"){
    if (methods[2] == TRUE){
      output$Method.2.data = Tm_data
      output$Method.2.fit <- Tm_vs_lnCt_fit
    }
  }
  if (methods[3] == TRUE){
    output$Method.3.data = gfit_data
    output$Method.3.fit <- gfit
  }
  output <- output
}
