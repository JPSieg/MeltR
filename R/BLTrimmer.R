#'Auto-trim absorbance melting curves to obtain more accurate meltR.A fits
#'
#'Automates #'
#'
#'@param meltR.A.fit data_frame containing absorbance melting data
#'@param n.combinations The blank sample for background subtraction, or a list of blanks to apply to different samples for background subtraction. "none" to turn off background subtraction. If there is a single blank in the data set, The identity of the blank, for example, blank = 1 or blanke = "blank". If there are multiple blanks in the data, blank = list(c("blank 1", "Sample 1"), c("blank 2", "sample 2")) and so on. Sample identifiers should be what they are in the data frame. If you need to figure out what the sample identifiers are, run unique(df$Sample), where df is the name of the R data frame you are using, in your R console.
#'@param n.ranges A vector containing the Nucleic acid type and the sequences you are fitting.
#'@param range.step The molecular model you want to fit. Options: "Monomolecular.2State", "Monomolecular.3State", "Heteroduplex.2State", "Homoduplex.2State".
#'@param no.trim.range The thermodynamic model you want to fit. Options: "VantHoff". Default = "VantHoff".
#'@param parallel The temperature used to calculate the NucAcid concentration. Default = 90.
#'@param n.core Option to only fit certain temperature ranges for melting curves. Either a vector or a list. If this is set to a vector, meltR.A will only fit temperatures in this range for all melting curves Example = c(17, 75). If set to a list of vectors, meltR.A will change what values are fit for each curve. Example, list(c(0,100), c(17,75), .... , c(0,100)). The length of this list has to be the equal to the number of samples that will be fit.
#'@param Save_results what methods do you want to use to fit data. Default = c(TRUE, TRUE, TRUE). Can be true or false. Note, method 1 must be set to TRUE or the subsequent steps will not work.
#'@param file_path What results to save. Options: "all" to save PDF plots and ".csv" formated tables of parameters, "some" to save ".csv" formated tables of parameters, or "none" to save nothing.
#'@param file_prefix Prefix that you want on the saved files.
#'@return A list of data frames containing parameters from the fits and data for plotting the results with ggplot2.
#' @export
BLTrimmer = function(meltR.A.fit,
                   n.combinations = 100,
                   n.ranges = 5,
                   range.step = 5,
                   no.trim.range = c(0.2, 0.8),
                   parallel = "none",
                   n.core = 1,
                   Save_results = "none",
                   file_path = getwd(),
                   file_prefix = "BLTrimmer"){
  ####Determine Mmodel and Tmodel####

  Tmodel = meltR.A.fit$BLTrimmer.data[[3]]
  Mmodel = meltR.A.fit$BLTrimmer.data[[2]]

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
      f = function(H, S, Temperature){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        model <- Mmodels$Monomolecular.2State(K = K)
        return(model)
      }
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
      f = function(H, S, Temperature, Ct){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        model <- Mmodels$Heteroduplex.2State(K = K, Ct = Ct)
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
      f = function(H, S, Temperature, Ct){
        K <- G_VH(H = H, S = S, Temperature = Temperature)
        model <- Mmodels$Homoduplex.2State(K = K, Ct = Ct)
        return(model)
      }
      calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(1/Ct)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))}
      calcG = function(H, Tm, Ct){ H - (310.15*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct)))) }
      calcG.SE = function(H, Tm, Ct, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct))))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2)}
    }
  }

  ####Get Cts for each sample####

  Samples = unique(meltR.A.fit$Method.1.data$Sample)

  list.df.raw = {}

  for (i in 1:length(Samples)){
    list.df.raw[[i]] = subset(meltR.A.fit$BLTrimmer.data[[1]], Sample == Samples[i])
    list.df.raw[[i]]$Ct = exp(meltR.A.fit$Method.2.data$lnCt[i])
  }

  ####Calculate the fraction unfolded for each reading####

  for (i in 1:length(Samples)){
    list.df.raw[[i]]$f = f(coef(meltR.A.fit$Method.3.fit)[1],
                           coef(meltR.A.fit$Method.3.fit)[2],
                           list.df.raw[[i]]$Temperature,
                           list.df.raw[[i]]$Ct)
  }

  ####Find the no trim range####

  list.no.trim = {}

  for (i in 1:length(Samples)){
    df = list.df.raw[[i]]
    high = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[1])])]
    low = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[2])])]
    list.no.trim[[i]] = c(low, high)
  }

  ####Create a list of baseline ranges for each sample####

  list.df.ranges = {}

  for (i in 1:length(Samples)){
    df = list.df.raw[[i]]
    v.samples = c()
    lows = c()
    highs = c()
    Windows = c()
    for (j in 1:n.ranges){
      v.samples[j] = df$Sample[1]
      lows[j] =  list.no.trim[[i]][1] - range.step*j
      highs[j] =  list.no.trim[[i]][2] + range.step*j
      Windows[j] = paste(v.samples[j], lows[j], "to", highs[j], sep = "#")
    }
    list.df.ranges[[i]] = Windows
  }

  #####Combine all combinations####

  baselines = expand.grid(list.df.ranges)
  baselines = baselines[runif(n.combinations, min = 1, max = nrow(baselines)),]

  ####Print number of baselines ####

  print(paste("You are trying to test", nrow(baselines), "baseline combinations"))
  print("Do you think this is possible?")

  ####Get starting parameters####

  list.fit.start.indv = {}

  for (j in 1:length(Samples)){
    list.fit.start.indv[[j]] <- list(H = coef(meltR.A.fit$Method.1.fit[[j]])[1],
                                     Tm = coef(meltR.A.fit$Method.1.fit[[j]])[2],
                                     mED = coef(meltR.A.fit$Method.1.fit[[j]])[3],
                                     bED = coef(meltR.A.fit$Method.1.fit[[j]])[4],
                                     mESS = coef(meltR.A.fit$Method.1.fit[[j]])[5],
                                     bESS = coef(meltR.A.fit$Method.1.fit[[j]])[6])
  }


  ####Function to fit baselines####

  BL.fitter = function(BL.df){

    if (parallel == "none"){
      pb = txtProgressBar(min = 1, max = nrow(baselines), initial = 1, style = 3)
    }

    #Define variables
    list.fit.start = {}
    frac.dH.error = c()
    dH1 = c()
    dS1 = c()
    dH2 = c()
    dS2 = c()

    for (i in 1:nrow(baselines)){
      #print(i)
      if (parallel == "none"){
        setTxtProgressBar(pb, i)
      }
      tryCatch({
        list.df = {}

        for (j in 1:length(Samples)) {
          #print(j)
          window = BL.df[i,j]
          window.vector = strsplit(as.character(window), "#")[[1]]
          low = as.numeric(window.vector[2])
          high = as.numeric(window.vector[4])
          df = subset(list.df.raw[[j]], Temperature <= high)
          df = subset(df, Temperature >= low)
          df$Sample = j
          list.df[[j]] = df
        }

        ####Method 1#####

        a <-{}
        fit <- {}
        indvfits.H <- c()
        indvfits.S <- c()
        Ind.model = Model
        mED = c()
        bED = c()
        mESS = c()
        bESS = c()
        for (j in c(1:length(Samples))){
          a[[j]] <- list.df[[j]]
          if (Mmodel == "Monomolecular.2State"){
            fit[[j]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature),
                            data = a[[j]],
                            start = list.fit.start.indv[[j]],
                            trace = FALSE,
                            nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
          }
          if (Mmodel == "Heteroduplex.2State"){
            fit[[j]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature, Ct),
                            data = a[[j]],
                            start = list.fit.start.indv[[j]],
                            trace = FALSE,
                            nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
          }
          if (Mmodel == "Homoduplex.2State"){
            fit[[j]] <- nls(Absorbance ~ Model(H, Tm, mED, bED, mESS, bESS, Temperature, Ct),
                            data = a[[j]],
                            start = list.fit.start.indv[[j]],
                            trace = FALSE,
                            nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))
          }
          indvfits.H[j] = coef(fit[[j]])[1]
          indvfits.S[j] = calcS(coef(fit[[j]])[1], coef(fit[[j]])[2], a[[j]]$Ct[1])
          mED[j] = coef(fit[[j]])[3]
          bED[j] = coef(fit[[j]])[4]
          mESS[j] = coef(fit[[j]])[5]
          bESS[j] = coef(fit[[j]])[6]
        }

        list.fit.start[[i]] = list.fit.start.indv

        dH1[i] = mean(indvfits.H)
        dS1[i] = mean(indvfits.S)

        ####Method 2####

        Tm_range <- {}
        Tm_fit <- {}
        Tm <- c()
        lnCt <- c()
        for (j in c(1:length(Samples))){
          #calculate f at each temp
          a[[j]]$f = (a[[j]]$Absorbance - (mESS[j]*a[[j]]$Temperature + bESS[j]))/((mED[j]*a[[j]]$Temperature + bED[j]) - (mESS[j]*a[[j]]$Temperature + bESS[j]))
          Tm_range[[j]] <- data.frame( #Pull out data where f >= 0.4 and f <= 0.6
            "Temperature" = a[[j]]$Temperature[which(a[[j]]$f >= 0.4)][which(a[[j]]$f[which(a[[j]]$f >= 0.4)] <= 0.6)],
            "f" = a[[j]]$f[which(a[[j]]$f >= 0.4)][which(a[[j]]$f[which(a[[j]]$f >= 0.4)] <= 0.6)]
          )
          Tm_fit[[j]] <- lm(f~Temperature, data = Tm_range[[j]]) #fit data where f >= 0.4 and f <= 0.6
          Tm[j] <- (0.5 - coef(Tm_fit[[j]])[1])/coef(Tm_fit[[j]])[2]
          lnCt[j] <- log(a[[j]]$Ct[1])
        }
        Tm_data <- data.frame("lnCt" = lnCt, "Tm" = Tm)
        Tm_data$invT <- 1/(Tm_data$Tm + 273.15)
        Tm_vs_lnCt_fit <- nls(invT ~ TmModel(H, S, lnCt),
                              data = Tm_data,
                              start = list(H = dH1[i], S = dS1[i]))

        dH2[i] = coef(Tm_vs_lnCt_fit)[1]
        dS2[i] = coef(Tm_vs_lnCt_fit)[2]

        frac.dH.error[i] = abs(dH1[i] - dH2[i])/abs(mean(dH1[i], dH2[i]))

      },
      error = function(e){})
    }
    BL.df$dH1 = dH1
    BL.df$dH2 = dH2
    BL.df$dS1 = dS1
    BL.df$dS2 = dS2
    BL.df$frac.dH.error = frac.dH.error
    output = BL.df
  }

  ####Unparallelized####

    print(paste("Fitting", nrow(baselines), "combinations of", n.ranges, "different baselines per sample"))

    if (parallel == "none"){
      stime <- system.time({

       baselines = BL.fitter(baselines)

      })
    }

    ####Paralleled####

    list.result = {}

    if (parallel == "on"){
      stime <- system.time({

        list.BL.df = split(baselines, rep(1:n.core, length.out = nrow(baselines), each = ceiling(nrow(baselines)/n.core)))

        c1 = parallel::makeCluster(2)
        doParallel::registerDoParallel(cores = n.core)

        list.result = foreach::foreach(k = 1:length(list.BL.df))  %dopar% {
          list.BL.df[[k]] = BL.fitter(list.BL.df[[k]])
        }

        baselines = list.result[[1]]

        if (length(list.BL.df) > 1){
          for (k in 2:length(list.BL.df)){
            baselines = rbind(baselines, list.result[[k]])
          }
        }

      })
    }

    print(paste("Testing baselines took", stime[3], "seconds"))

  ####Grab best baseline####

  df.raw = meltR.A.fit$BLTrimmer.data[[1]]

  df.opt = baselines[which.min(baselines$frac.dH.error),]

  fitstart = list.fit.start.indv

  v.sample = c()

  list.T.range = {}
  list.df.trim = {}
  list.df.no.trim = {}

  for (i in 1:length(Samples)){
    window = df.opt[,i]
    window.vector = strsplit(as.character(window), "#")[[1]]
    low = as.numeric(window.vector[2])
    high = as.numeric(window.vector[4])
    list.T.range[[i]] = c(low, high)
    list.df.trim[[i]] = subset(df.raw, Sample == Samples[i])
    list.df.trim[[i]] = subset(list.df.trim[[i]], Temperature >= low)
    list.df.trim[[i]] = subset(list.df.trim[[i]], Temperature <= high)
    list.df.no.trim[[i]] = subset(df.raw, Sample == Samples[i])
    list.df.no.trim[[i]] = subset(list.df.no.trim[[i]], Temperature >= list.no.trim[[i]][1])
    list.df.no.trim[[i]] = subset(list.df.no.trim[[i]], Temperature <= list.no.trim[[i]][2])
  }

  df.trim = list.df.trim[[1]]
  df.no.trim = list.df.no.trim[[1]]

  for (i in 2:length(Samples)){
    df.trim = rbind(df.trim, list.df.trim[[i]])
    df.no.trim = rbind(df.no.trim, list.df.no.trim[[i]])
  }

  ####Make plots####

  if (Save_results == "all"){

    #Histogram

    pdf(paste(file_path, "/", file_prefix, "_histogram.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    hist(baselines$frac.dH.error,
         xlab = "|dH1 - dH2|/Average dH",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8)
    dev.off()

    #dH vrs. frac error

    pdf(paste(file_path, "/", file_prefix, "_dH_vs_frac_error.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(baselines$dH1 ~ log10(baselines$frac.dH.error),
         xlab = "|dH1 - dH2|/Average dH",
         ylab = "dH1 (black) or dH2 (red)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                  max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)))
    par(new=T)
    plot(baselines$dH2 ~ log10(baselines$frac.dH.error), col = "red",
         xlab = "|dH1 - dH2|/Average dH",
         ylab = "dH1 (black) or dH2 (red)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                  max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)))
    dev.off()

    #Trimmed baseline plot

    pdf(paste(file_path, "/", file_prefix, "_trimmed_baseline_plot.pdf", sep = ""),
        width = 3, height = 3, pointsize = 0.25)
    plot(df.raw$Absorbance ~ df.raw$Temperature,
         xlab = Temperature ~ (degree ~ C), ylab = "Absorbance",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(df.raw$Absorbance),
                  max(df.raw$Absorbance)),
         xlim = c(min(df.raw$Temperature),
                  max(df.raw$Temperature)))
    par(new=T)
    plot(df.trim$Absorbance ~ df.trim$Temperature, col = "red",
         xlab = "", ylab = "",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(df.raw$Absorbance),
                  max(df.raw$Absorbance)),
         xlim = c(min(df.raw$Temperature),
                  max(df.raw$Temperature)))
    par(new=T)
    plot(df.no.trim$Absorbance ~ df.no.trim$Temperature, col = "blue",
         xlab = "", ylab = "",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(df.raw$Absorbance),
                  max(df.raw$Absorbance)),
         xlim = c(min(df.raw$Temperature),
                  max(df.raw$Temperature)))
    dev.off()

  }

  ####Fit auto trimmed data with meltR.A####

  print("Using autotrimmed baselines in MeltR")

  meltR.A.fit = meltR.A(meltR.A.fit$meltR.A.settings[[1]],
                        meltR.A.fit$meltR.A.settings[[2]],
                        meltR.A.fit$meltR.A.settings[[3]],
                        meltR.A.fit$meltR.A.settings[[4]],
                        list.T.range,
                        meltR.A.fit$meltR.A.settings[[6]],
                        meltR.A.fit$meltR.A.settings[[7]],
                        meltR.A.fit$meltR.A.settings[[8]],
                        Save_results,
                        file_prefix,
                        file_path,
                        fitstart)

  output = list(list.T.range,
                meltR.A.fit,
                stime)
  names(output) = c("Optimum.T.ranges", "Optimum.fit", "System.time")

  output = output

}
