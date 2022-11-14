#'Auto-trim absorbance melting curves to obtain more accurate meltR.A fits
#'
#'Automates absorbance data baseline trimming and testing. Generates random permutations of trimmed baseline combinations, fits the data using Method 1 and Method 2, compares the enthalpy difference, and chooses an optimum set of trimmed absorbance melting curves.
#'
#'@param meltR.A.fit A object produced by fitting data with MeltR.A.
#'@param Trim.method Method for trimming baselines. "fixed" to use the same baseline lengths for each curve or "floating" to use different baseline lengths for each curve. Default = "floating".
#'@param Assess.method Method for assessing fits from trimmed baseline. Options are integers 1, 2, and 3. 1 maximizes agreement between individual fits. 2 maximizes agreement between the average of individual fits and the average of the 1/Tm versus lnCt method. 3 takes both 1 and 2 into account. Default = 3.
#'@param n.combinations Number of baseline combinations to test using the float method. The program will produce n.ranges^Samples combinations of baselines. It will require a large amount of computational time to test these. In general, testing 1000 combinations will produce a reliable result (plus or minus 5% in terms on enthalpy). For an exhaustive testing, set this parameter to n.ranges^Samples.
#'@param n.ranges.float Number of trimmed baselines to generate per sample using the float method. It is not recommended to increase this parameter past 6 because of how long it will take the computer to generate all of these combinations.
#'@param range.step.float Temperature difference between each range that is generated for each sample using the float method. Default = 5 deg Celsius works well.
#'@param n.ranges.fixed Number of baseline ranges for the fixed method.
#'@param range.step.fixed. Temperature difference between each range that is generated using the fixed method.
#'@param no.trim.range Determines the range where the absorbance data will not be trimmed. By default, no.trim.range = c(0.2, 0.8), meaning that the data will not be trimmed at a mode fraction double stranded greater than 0.2 and less than 0.8. Determined using the global fit from the input object created by meltR.A.
#'@param quantile.threshold Threshold for assessing the best baseline combinations
#'@param parallel This is argument is not exposed to end users. Contact Jacob Sieg if you are interested. to "on" if you want to use multiple cores to fit baselines. You will need to run library(doParallel) to load doParallel and its dependencies. You will also need to specify the number of cores in the n.cores argument, which should not exceed the number of cores your computer has. There will be no benefit for running in parallel mode for a single core machine. Default = "none".
#'@param n.core Associated with n.core and not exposed to end users. Contact Jacob Sieg if you are interested. How many cores do you want to designate for this task.
#'@param Save_results "all" to save results to the disk or "none" to not save results to the disk.
#'@param file_path A path to the folder you want your results saved in.
#'@param file_prefix Prefix that you want on the saved files.
#'@param memory.light Keeps the BLTrimmer from blowing up memory usage. If TRUE, only passes the objects that the BLtrimmer needs to run to the memory. If FALSE passes a much more extensive number of statistics to the memory.
#'@param Silent Set to TRUE to not print results
#'@return  A BLTrimmer fit object containing a list of data objects for advanced analysis. Set memory.light = FALSE to get the full object.
#' \itemize{
#'   \item 1. Baseline.data - A data frame containing the temperature range for each sample with the dH and dS for methods 1 and 2. Also contains error between each method. This object is not returned when memory.light = TRUE.
#'   \item 2. List.T.ranges - A list containing the temperature range for each sample and baseline combination tested. This object is not returned when memory.light = TRUE.
#'   \item 3. List.fits - A list containing the meltR.A fit object for each baseline combination that was passed to the optimum baseline ensemble. This object is not returned when memory.light = TRUE.
#'   \item 4. Fit.summaries - A data frame containing thermodynamic parameters from a meltR.A fit object for each baseline combination that was passed to the optimum baseline ensemble. This object is not returned when memory.light = TRUE.
#'   \item 5. Ensemble.energies - A data frame containing the thermodynamic parameters from an analysis of the optimum ensemble energies. This object is returned when memory.light = TRUE.
#'   \item 6. Fractional.error.between.methods - A data frame containing percent error between methods from an analysis of the optimum ensemble energies. This object is returned when memory.light = TRUE.
#'   \item 7. System.time - A list containing the system time used by the BLTrimmer run. The first element is how long the fast analysis of many baseline combinations took. The second element is how long the fitting of the optimum baseline ensemble by meltR.A took. This object is returned when memory.light = TRUE.
#'   }
#' @export
BLTrimmer = function(meltR.A.fit,
                     Trim.method = "floating",
                     Assess.method = 3,
                     n.combinations = 1000,
                     n.ranges.float = 5,
                     range.step.float = 5,
                     n.ranges.fixed = 40,
                     range.step.fixed = 0.5,
                     no.trim.range = c(0.1, 0.9),
                     quantile.threshold = 0.25,
                     Save_results = "none",
                     file_path = getwd(),
                     file_prefix = "BLTrimmer",
                     memory.light = TRUE,
                     Silent = FALSE){
  ####Potential arguments for parallel computing on multicore machines####

  parallel = "none" #This will keep the parallel code from running on accident
  n.core = 1

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
    list.df.raw[[i]]$Ct = meltR.A.fit$Method.1.indvfits$Ct[i]
  }

  ####Calculate the fraction unfolded for each reading####

  for (i in 1:length(Samples)){
    if (Mmodel == "Monomolecular.2State"){
      list.df.raw[[i]]$f = f(coef(meltR.A.fit$Method.3.fit)[1],
                             coef(meltR.A.fit$Method.3.fit)[2],
                             list.df.raw[[i]]$Temperature)
    }else{
      list.df.raw[[i]]$f = f(coef(meltR.A.fit$Method.3.fit)[1],
                             coef(meltR.A.fit$Method.3.fit)[2],
                             list.df.raw[[i]]$Temperature,
                             list.df.raw[[i]]$Ct)
    }
  }

  ####Find the no trim range####

  list.no.trim = {}

  for (i in 1:length(Samples)){
    df = list.df.raw[[i]]
    high = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[1])])]
    low = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[2])])]
    if (length(low) == 0){
      low = df$Temperature[which.min(df$Temperature)]
    }
    if (length(high) == 0){
      high = df$Temperature[which.max(df$Temperature)]
    }

    list.no.trim[[i]] = c(low, high)
  }

  ####Create a list of baseline ranges for each sample####

  ####Trim method 1####

  if (Trim.method == "fixed"){
    list.df.ranges = {}

    for (i in 1:length(Samples)){
      df = list.df.raw[[i]]
      v.samples = c()
      lows = c()
      highs = c()
      Windows = c()
      for (j in 1:n.ranges.fixed){
        v.samples[j] = df$Sample[1]
        lows[j] =  list.no.trim[[i]][1] - range.step.fixed*j
        highs[j] =  list.no.trim[[i]][2] + range.step.fixed*j
        Windows[j] = paste(v.samples[j], lows[j], "to", highs[j], sep = "#")
      }
      list.df.ranges[[i]] = Windows
    }
    baselines = matrix(nrow = n.ranges.fixed, ncol = length(Samples))

    for (i in 1:n.ranges.fixed){
      #print(i)
      for (j in 1:length(Samples)){
        #print(j)
        baselines[i,j] = list.df.ranges[[j]][i]
      }
    }
  }

  ####Trim method 2####

  if (Trim.method == "floating"){
    list.df.ranges = {}

    for (i in 1:length(Samples)){
      df = list.df.raw[[i]]
      v.samples = c()
      lows = c()
      highs = c()
      Windows = c()
      for (j in 1:n.ranges.float){
        v.samples[j] = df$Sample[1]
        lows[j] =  list.no.trim[[i]][1] - range.step.float*j
        highs[j] =  list.no.trim[[i]][2] + range.step.float*j
        Windows[j] = paste(v.samples[j], lows[j], "to", highs[j], sep = "#")
      }
      list.df.ranges[[i]] = Windows
    }

    baselines = expand.grid(list.df.ranges)
    baselines = baselines[runif(n.combinations, min = 1, max = nrow(baselines)),]

  }
  ####Print number of baselines ####

  if(Silent){}else{
    print(paste("You are trying to test", nrow(baselines), "baseline combinations"))
    print("Do you think this is possible?")
  }

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
      if (Silent){}else{pb = txtProgressBar(min = 1, max = nrow(baselines), initial = 1, style = 3)}
    }

    #Define variables
    list.fit.start = {}
    frac.dH1.error = c()
    frac.dH1.dH2.error = c()
    dH1 = c()
    dS1 = c()
    dH2 = c()
    dS2 = c()


    for (i in 1:nrow(baselines)){
      #print(i)
      if (parallel == "none"){
        if (Silent){}else{ setTxtProgressBar(pb, i)}
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
        indvfits.Tm = c()
        indvfits.Ct = c()
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

          if (Mmodel == "Monomolecular.2State"){
            indvfits.S[j] = calcS(coef(fit[[j]])[1], coef(fit[[j]])[2])
          }else{
            indvfits.S[j] = calcS(coef(fit[[j]])[1], coef(fit[[j]])[2], a[[j]]$Ct[1])
          }

          indvfits.Tm[j] = coef(fit[[j]])[2]
          indvfits.Ct[j] = a[[j]]$Ct[1]
          mED[j] = coef(fit[[j]])[3]
          bED[j] = coef(fit[[j]])[4]
          mESS[j] = coef(fit[[j]])[5]
          bESS[j] = coef(fit[[j]])[6]
        }

        list.fit.start[[i]] = list.fit.start.indv

        dH1[i] = mean(indvfits.H)
        dS1[i] = mean(indvfits.S)

        ####Method 2####

        if (Mmodel == "Monomolecular.2State"){
          dH2[i] = NA
          dS2[i] = NA
          frac.dH1.dH2.error[i] = NA
          frac.dH1.error[i] = abs(sd(indvfits.H)/mean(indvfits.H))
        }else{

          lnCt = log(indvfits.Ct)
          Tm = indvfits.Tm
          Tm_data = data.frame(lnCt, Tm)
          Tm_data$invT <- 1/(Tm_data$Tm + 273.15)
          Tm_vs_lnCt_fit <- nls(invT ~ TmModel(H, S, lnCt),
                                data = Tm_data,
                                start = list(H = dH1[i], S = dS1[i]))

          dH2[i] = coef(Tm_vs_lnCt_fit)[1]
          dS2[i] = coef(Tm_vs_lnCt_fit)[2]

          frac.dH1.dH2.error[i] = abs(dH1[i] - dH2[i])/abs(mean(dH1[i], dH2[i]))
          frac.dH1.error[i] = abs(sd(indvfits.H)/mean(indvfits.H))
        }


      },
      error = function(e){})
    }
    BL.df=data.frame(BL.df)
    BL.df$dH1 = dH1
    BL.df$dH2 = dH2
    BL.df$dS1 = dS1
    BL.df$dS2 = dS2
    BL.df$frac.dH1.error = frac.dH1.error
    BL.df$frac.dH1.dH2.error =  frac.dH1.dH2.error
    output = BL.df
  }

  ####Unparallelized####

  if (Silent){}else{
    if (Trim.method == "fixed"){
      print(paste("Fitting", nrow(baselines), "baseline combinations"))
    }
    if (Trim.method == "floating"){
      print(paste("Fitting", nrow(baselines), "combinations of", n.ranges.float, "different baselines per sample"))
    }
  }

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

  ####Quantilize data####

  quantiles.dH1.error = quantile(baselines$frac.dH1.error, seq(0, 1, length.out = nrow(baselines)))

  if (Mmodel == "Monomolecular.2State"){
    quantiles.dH1.dH2.error = NA
  }else{
    quantiles.dH1.dH2.error = quantile(baselines$frac.dH1.dH2.error, seq(0, 1, length.out = nrow(baselines)))
  }


  dH1.error.quantile = c()
  dH1.dH2.error.quantile = c()

  for (i in 1:nrow(baselines)) {
    dH1.error.quantile[i] = as.numeric(gsub("%", "", names(quantiles.dH1.error)[which.min(abs(quantiles.dH1.error - baselines$frac.dH1.error[i]))]))/100
    if (Mmodel == "Monomolecular.2State"){
      dH1.dH2.error.quantile[i] = NA
    }else{
      dH1.dH2.error.quantile[i] = as.numeric(gsub("%", "", names(quantiles.dH1.error)[which.min(abs(quantiles.dH1.dH2.error - baselines$frac.dH1.dH2.error[i]))]))/100
    }

  }

  if (Mmodel == "Monomolecular.2State"){
    error.distance = NA
  }else{
    error.distance = sqrt(dH1.dH2.error.quantile^2 + dH1.error.quantile^2)
  }

  baselines$dH1.error.quantile = dH1.error.quantile
  baselines$dH1.dH2.error.quantile = dH1.dH2.error.quantile
  baselines$error.distance =  error.distance

  ####Pull out top 10% of best values based on some criterion####

  if (Mmodel == "Monomolecular.2State"){
    Assess.method = 1
  }

  #Method 1 agrees with method 1

  n.best = ceiling(nrow(baselines)*quantile.threshold)

  baselines = baselines[order(baselines$dH1.error.quantile),]

  df.m1 = baselines[1:n.best,]

  #Method 2 agrees with method 1

  n.best = ceiling(nrow(baselines)*quantile.threshold)

  baselines = baselines[order(baselines$dH1.dH2.error.quantile),]

  df.m2 = baselines[1:n.best,]

  #Method 1 agrees with method 1 and Method 2 agrees with method 1

  n.best = ceiling(nrow(baselines)*quantile.threshold)

  baselines = baselines[order(baselines$error.distance),]

  df.m1.M2 = baselines[1:n.best,]

  if (Assess.method == 1){
    df.best = df.m1
    color = "orange"
  }
  if (Assess.method == 2){
    df.best = df.m2
    color = "cyan"
  }
  if (Assess.method == 3){
    df.best = df.m1.M2
    color = "blue"
  }

  ####Make plots####

  if (Save_results == "all"){

    pdf(paste(file_path, "/", file_prefix, "_baseline_ensemble_analysis.pdf", sep = ""),
        width = 9, height = 6, pointsize = 0.25)


    layout_matrix_1 <- matrix(1:6, ncol = 3) # Define position matrix
    layout(layout_matrix_1)

    #Histogram 1 agree 1

    hist(baselines$frac.dH1.error,
         xlab = "SD(dH1)/Average dH",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         main = "Standard deviation of Method 1 dH")
    abline(v = df.best$frac.dH1.error, col = adjustcolor(color, alpha = 0.3))

    #dH1 error vrs. dH1 error

    plot(baselines$dH1 ~ baselines$frac.dH1.error,
         xlab = "SD(dH1)/Average dH",
         ylab = "dH1 (black)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
         ylim = c(min(c(baselines$dH1), na.rm = TRUE),
                  max(c(baselines$dH1), na.rm = TRUE)),
         xlim = c(min(c(baselines$frac.dH1.error), na.rm = TRUE),
                  max(c(baselines$frac.dH1.error), na.rm = TRUE)),
         main = "Standard deviation of Method 1 dH")
    par(new=T)
    plot(df.best$dH1 ~ df.best$frac.dH1.error,
         xlab = "SD(dH1)/Average dH",
         ylab = "dH1 (black)",
         cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = adjustcolor(color),
         ylim = c(min(c(baselines$dH1), na.rm = TRUE),
                  max(c(baselines$dH1), na.rm = TRUE)),
         xlim = c(min(c(baselines$frac.dH1.error), na.rm = TRUE),
                  max(c(baselines$frac.dH1.error), na.rm = TRUE)))

    #Histogram 1 agree 2

    if (Mmodel == "Monomolecular.2State"){
      hist(0,
           xlab = "|dH1 - dH2|/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           main = "Method 1 dH agrees with Method 2 dH")
    }else{
      hist(baselines$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           main = "Method 1 dH agrees with Method 2 dH")
      abline(v = df.best$frac.dH1.dH2.error, col = adjustcolor(color, alpha = 0.3))
    }


    #dH vrs. dH1 to dH2 error

    if (Mmodel == "Monomolecular.2State"){
      plot(0 ~ 1,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "dH1 (black) or dH2 (red)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           main = "Method 1 dH agrees with Method 2 dH")
    }else{
      plot(baselines$dH1 ~ baselines$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "dH1 (black) or dH2 (red)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                    max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)),
           main = "Method 1 dH agrees with Method 2 dH")
      par(new=T)
      plot(baselines$dH2 ~ baselines$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "dH1 (black) or dH2 (red)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = "red",
           ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                    max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)))
      par(new=T)
      plot(df.best$dH1 ~ df.best$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "dH1 (black) or dH2 (red)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = adjustcolor(color),
           ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                    max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)))
      par(new=T)
      plot(df.best$dH2 ~ df.best$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "dH1 (black) or dH2 (red)",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = adjustcolor(color),
           ylim = c(min(c(baselines$dH1, baselines$dH2), na.rm = TRUE),
                    max(c(baselines$dH1, baselines$dH2), na.rm = TRUE)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)))
    }



    #dH1 to dH2 error quantile verseus dH1 error quantile

    if (Mmodel == "Monomolecular.2State"){
      plot(0 ~ 1,
           xlab = "Quantile |dH1 - dH2|/Average dH",
           ylab = "Quantile SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           main = "Error distance")
    }else{
      plot(baselines$dH1.error.quantile ~ baselines$dH1.dH2.error.quantile,
           xlab = "Quantile |dH1 - dH2|/Average dH",
           ylab = "Quantile SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           ylim = c(0,1),
           xlim = c(0,1),
           main = "Error distance")
      par(new=T)
      plot(df.best$dH1.error.quantile ~ df.best$dH1.dH2.error.quantile,
           xlab = "Quantile |dH1 - dH2|/Average dH",
           ylab = "Quantile SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = adjustcolor(color),
           ylim = c(0,1),
           xlim = c(0,1))
      for (i in 1:nrow(df.best)){
        segments(x0 = 0,
                 y0 = 0,
                 x1 = df.best$dH1.dH2.error.quantile[i],
                 y1 = df.best$dH1.error.quantile[i],
                 col = adjustcolor(color, alpha = 0.3))
      }
    }



    #dH1 error verseus dH1 to dH2 error

    if (Mmodel == "Monomolecular.2State"){
      plot(0 ~ 1,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           main = "Error distance")
    }else{
      plot(baselines$frac.dH1.error ~ baselines$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8,
           ylim = c(min(baselines$frac.dH1.error),
                    max(baselines$frac.dH1.error)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)),
           main = "Error distance")
      par(new=T)
      plot(df.best$frac.dH1.error ~ df.best$frac.dH1.dH2.error,
           xlab = "|dH1 - dH2|/Average dH",
           ylab = "SD(dH1)/Average dH",
           cex.lab = 1.5, cex.axis = 1.25, cex = 0.8, col = adjustcolor(color),
           ylim = c(min(baselines$frac.dH1.error),
                    max(baselines$frac.dH1.error)),
           xlim = c(min(baselines$frac.dH1.dH2.error),
                    max(baselines$frac.dH1.dH2.error)))
    }

    dev.off()

  }

  ####Set up meltR.A####

  if (Silent){}else{
    print("")
    print("Using autotrimmed baselines in meltR.A")
  }

  list.meltR.A.fit = {}
  list.df.results = {}
  list.baselines = {}

  if (Silent){}else{pb = txtProgressBar(min = 1, max = nrow(df.best), initial = 1, style = 3)}

  stime2 <- system.time({

  for (k in 1:nrow(df.best)){

    if(Silent){}else{setTxtProgressBar(pb, k)}

    df.raw = meltR.A.fit$BLTrimmer.data[[1]]

    df.opt = df.best[k,]

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

    meltR.A.fit = meltR.A(data_frame = meltR.A.fit$meltR.A.settings[[1]],
                          blank = meltR.A.fit$meltR.A.settings[[2]],
                          NucAcid = meltR.A.fit$meltR.A.settings[[3]],
                          concT = meltR.A.fit$meltR.A.settings[[4]],
                          fitTs = list.T.range,
                          methods = meltR.A.fit$meltR.A.settings[[6]],
                          Mmodel = meltR.A.fit$meltR.A.settings[[7]],
                          Tmodel = meltR.A.fit$meltR.A.settings[[8]],
                          auto.trimmed = fitstart,
                          Silent = TRUE)
    if (memory.light){
      list.meltR.A.fit[[k]] = NA
      list.df.results[[k]] = meltR.A.fit$Summary
      list.baselines[[k]] = NA
    }else{
      list.meltR.A.fit[[k]] = meltR.A.fit
      list.df.results[[k]] = meltR.A.fit$Summary
      list.baselines[[k]] = list.T.range
    }
  }
  })

  df.result = list.df.results[[1]]

  for (i in 2:length(list.meltR.A.fit)){
    df.result = rbind(df.result, list.df.results[[i]])
  }

  df1.data = subset(df.result, Method == "1 individual fits")

  df2.data = subset(df.result, Method == "2 Tm versus ln[Ct]")

  df3.data = subset(df.result, Method == "3 Global fit")

  if (Trim.method == "fixed"){
    df.1 = data.frame("Method" = "1 individual fits",
                      "dH" = round(mean(df1.data$dH), digits = 2),
                      "SE.dH" = round(sd(df1.data$dH), digits = 2),
                      "dS" = round(mean(df1.data$dS), digits = 2),
                      "SE.dS"= round(sd(df1.data$dS), digits = 2),
                      "dG" = round(mean(df1.data$dG), digits = 2),
                      "SE.dG" = round(sd(df1.data$dG), digits = 2),
                      "Tm_at_0.1mM" = round(mean(df1.data$Tm_at_0.1mM), digits = 2),
                      "SE.Tm_at_0.1mM" = round(sd(df1.data$Tm_at_0.1mM), digits = 2))
    df.2 = data.frame("Method" = "2 Tm versus ln[Ct]",
                      "dH" = round(mean(df2.data$dH), digits = 2),
                      "SE.dH" = round(sd(df2.data$dH), digits = 2),
                      "dS" = round(mean(df2.data$dS), digits = 2),
                      "SE.dS"= round(sd(df2.data$dS), digits = 2),
                      "dG" = round(mean(df2.data$dG), digits = 2),
                      "SE.dG" = round(sd(df2.data$dG), digits = 2),
                      "Tm_at_0.1mM" = round(mean(df2.data$Tm_at_0.1mM), digits = 2),
                      "SE.Tm_at_0.1mM" = round(sd(df2.data$Tm_at_0.1mM), digits = 2))
    df.3 = data.frame("Method" = "3 Global fit",
                      "dH" = round(mean(df3.data$dH), digits = 2),
                      "SE.dH" = round(sd(df3.data$dH), digits = 2),
                      "dS" = round(mean(df3.data$dS), digits = 2),
                      "SE.dS"= round(sd(df3.data$dS), digits = 2),
                      "dG" = round(mean(df3.data$dG), digits = 2),
                      "SE.dG" = round(sd(df3.data$dG), digits = 2),
                      "Tm_at_0.1mM" = round(mean(df3.data$Tm_at_0.1mM), digits = 2),
                      "SE.Tm_at_0.1mM" = round(sd(df3.data$Tm_at_0.1mM), digits = 2))
  }else{
    df.1 = data.frame("Method" = "1 individual fits",
                      "dH" = round(mean(df1.data$dH), digits = 2),
                      "CI95.dH" = paste(round(quantile(df1.data$dH, 0.025), 2), "to", round(quantile(df1.data$dH, 0.975), 2)),
                      "dS" = round(mean(df1.data$dS), digits = 2),
                      "CI95.dS"= paste(round(quantile(df1.data$dS, 0.025), 2), "to", round(quantile(df1.data$dS, 0.975), 2)),
                      "dG" = round(mean(df1.data$dG), digits = 2),
                      "CI95.dG" = paste(round(quantile(df1.data$dG, 0.025), 2), "to", round(quantile(df1.data$dG, 0.975), 2)),
                      "Tm_at_0.1mM" = round(mean(df1.data$Tm_at_0.1mM), digits = 2),
                      "CI95.Tm_at_0.1mM" = paste(round(quantile(df1.data$Tm_at_0.1mM, 0.025), 2), "to", round(quantile(df1.data$Tm_at_0.1mM, 0.975), 2)))
    df.2 = data.frame("Method" = "2 Tm versus ln[Ct]",
                      "dH" = round(mean(df2.data$dH), digits = 2),
                      "CI95.dH" = paste(round(quantile(df2.data$dH, 0.025), 2), "to", round(quantile(df2.data$dH, 0.975), 2)),
                      "dS" = round(mean(df2.data$dS), digits = 2),
                      "CI95.dS"= paste(round(quantile(df2.data$dS, 0.025), 2), "to", round(quantile(df2.data$dS, 0.975), 2)),
                      "dG" = round(mean(df2.data$dG), digits = 2),
                      "CI95.dG" = paste(round(quantile(df2.data$dG, 0.025), 2), "to", round(quantile(df2.data$dG, 0.975), 2)),
                      "Tm_at_0.1mM" = round(mean(df2.data$Tm_at_0.1mM), digits = 2),
                      "CI95.Tm_at_0.1mM" = paste(round(quantile(df2.data$Tm_at_0.1mM, 0.025), 2), "to", round(quantile(df2.data$Tm_at_0.1mM, 0.975), 2)))
    df.3 = data.frame("Method" = "3 Global fit",
                      "dH" = round(mean(df3.data$dH), digits = 2),
                      "CI95.dH" = paste(round(quantile(df3.data$dH, 0.025), 2), "to", round(quantile(df3.data$dH, 0.975), 2)),
                      "dS" = round(mean(df3.data$dS), digits = 2),
                      "CI95.dS"= paste(round(quantile(df3.data$dS, 0.025), 2), "to", round(quantile(df3.data$dS, 0.975), 2)),
                      "dG" = round(mean(df3.data$dG), digits = 2),
                      "CI95.dG" = paste(round(quantile(df3.data$dG, 0.025), 2), "to", round(quantile(df3.data$dG, 0.975), 2)),
                      "Tm_at_0.1mM" = round(mean(df3.data$Tm_at_0.1mM), digits = 2),
                      "CI95.Tm_at_0.1mM" = paste(round(quantile(df3.data$Tm_at_0.1mM, 0.025), 2), "to", round(quantile(df3.data$Tm_at_0.1mM, 0.975), 2)))
  }


  df.final = rbind(df.1, df.2)
  df.final = rbind(df.final, df.3)

  df.error = data.frame("dH" = round(100*abs((min(df.final$dH) - max(df.final$dH))/mean(df.final$dH)), digits = 1),
                        "dS" = round(100*abs((min(df.final$dS) - max(df.final$dS))/mean(df.final$dS)), digits = 1),
                        "dG" = round(100*abs((min(df.final$dG) - max(df.final$dG))/mean(df.final$dG)), digits = 1))

  ####Assemble output####

  stime = list(stime, stime2)

  names(stime) = c("all", "meltR.A")

  if (memory.light){
    baselines = NA
    list.baselines = NA
    list.meltR.A.fit = NA
    df.result = NA
  }

  output = list(baselines,
                list.baselines,
                list.meltR.A.fit,
                df.result,
                df.final,
                df.error,
                stime)

  names(output) = c("Baseline.data","List.T.ranges", "List.fits", "Fit.summaries", "Ensemble.energies", "Fractional.error.between.methods", "System.time")

  if(Silent){}else{
    print("Ensemble energies")
    print(df.final)
    print("Fractional error between methods")
    print(df.error)
    print("dH and dG are in kcal/mol and dS is in cal/mol/K. Tms are in deg Celsius")
  }

  if (Save_results == "all"){
    write.csv(df.final, paste(file_path, "/", file_prefix, "_baseline_ensemble_analysis_results.csv", sep = ""),row.names = FALSE)
  }

  output = output

}
