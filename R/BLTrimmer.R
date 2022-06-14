library(MeltR)

df = df.abs.data

?df.abs.data
?meltR.A

list.T.range = list(c(15, 70), #Sample 2
                    c(15, 70), #Sample 3
                    c(15, 80), #Sample 4
                    c(15, 85), #Sample 5
                    c(20, 78), #Sample 6
                    c(20, 85), #Sample 7
                    c(20, 85), #Sample 8
                    c(25, 85), #Sample 9
                    c(25, 85)) #Sample 10

fit = meltR.A(df, blank = 1,
        NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
        Mmodel = "Heteroduplex.2State",
        fitTs = list.T.range)

raw.df = df
meltR.A.fit = fit
n.ranges = 5
range.step = 5
no.trim.range = c(0.2, 0.8)
parallel = "none"
n.core = 1

####Get Cts for each sample####

Samples = unique(meltR.A.fit$Method.1.data$Sample)

list.df.raw = {}

for (i in 1:length(Samples)){
  list.df.raw[[i]] = subset(raw.df, Sample == Samples[i])
  list.df.raw[[i]]$Ct = exp(meltR.A.fit$Method.2.data$lnCt[i])
}

####Calculate the fraction unfolded for each reading####

Samples = unique(meltR.A.fit$Method.3.data$Sample)

for (i in 1:length(Samples)){
  df = list.df.raw[[i]]
  fit = meltR.A.fit$Method.3.fit
  coeff = coef(fit)
  mED = coeff[which(names(coeff) == paste("mED", Samples[i], sep = ""))]
  bED = coeff[which(names(coeff) == paste("bED", Samples[i], sep = ""))]
  mESS = coeff[which(names(coeff) == paste("mESS", Samples[i], sep = ""))]
  bESS = coeff[which(names(coeff) == paste("bESS", Samples[i], sep = ""))]
  list.df.raw[[i]]$f = (df$Absorbance - (mESS*df$Temperature + bESS))/((mED*df$Temperature + bED) - (mESS*df$Temperature + bESS))
}

####Find the no trim range####

list.no.trim = {}

for (i in 1:length(Samples)){
  df = list.df.raw[[i]]
  high = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[1])])]
  low = df$Temperature[which.max(df$Temperature[which(df$f >= no.trim.range[2])])]
  list.no.trim[[i]] = c(low, high)
}

####Print number of baselines ####

print(paste("You are trying to test", n.ranges^length(Samples), "baseline combinations"))
print("Do you think this is possible?")

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

####Starting values for the global fit####

gfit_start = list(H = coeff[1],
                  S = coeff[2],
                  mED = coeff[3:(2+length(Samples))],
                  bED = coeff[(3+length(Samples)):((2+2*length(Samples)))],
                  mESS = coeff[((3+2*length(Samples))):((2+3*length(Samples)))],
                  bESS = coeff[((3+3*length(Samples))):((2+4*length(Samples)))])

####Unparallelized####

frac.dH.error = c()

print(paste("Fitting", n.ranges^length(Samples), "combinations of", n.ranges, "different baselines per sample"))

if (parallel == "none"){
  pb = txtProgressBar(min = 1, max = nrow(baselines), initial = 1, style = 3)
  for (i in 1:nrow(baselines)){
    #print(i)
    setTxtProgressBar(pb, i)
    tryCatch({
      list.df = {}

      for (j in 1:length(Samples)) {
        #print(j)
        window = baselines[i,j]
        window.vector = strsplit(as.character(window), "#")[[1]]
        low = as.numeric(window.vector[2])
        high = as.numeric(window.vector[4])
        df = subset(list.df.raw[[j]], Temperature <= high)
        df = subset(df, Temperature >= low)
        df$Sample = j
        list.df[[j]] = df
      }

      df = list.df[[1]]
      for (j in 2:length(Samples)) {
        df = rbind(df, list.df[[j]])
      }
      ####Method 3#####

      gfit <- nls(Absorbance ~ GModel(H, S, mED, bED, mESS, bESS, Sample, Temperature, Ct),
                  start = gfit_start,
                  data = df,
                  nls.control(tol = 5e-04, minFactor = 1e-10, maxiter = 50, warnOnly = TRUE))

      ####Method 2####

      coeff = coef(gfit)

      mED = coeff[3:(2+length(Samples))]
      bED = coeff[(3+length(Samples)):((2+2*length(Samples)))]
      mESS = coeff[((3+2*length(Samples))):((2+3*length(Samples)))]
      bESS = coeff[((3+3*length(Samples))):((2+4*length(Samples)))]

      a <- {}
      Tm_range <- {}
      Tm_fit <- {}
      Tm <- c()
      lnCt <- c()
      for (j in c(1:length(Samples))){
        a[[j]] <- subset(df, Sample == Samples[j]) #Pull out data
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
                            start = list(H = -70, S = -0.17))
      frac.dH.error[i] = abs(coef(Tm_vs_lnCt_fit)[1] - coef(gfit)[1])/abs(mean(coef(Tm_vs_lnCt_fit)[1], coef(gfit)[1]))
    },
             error = function(e){})
  }
}

hist(frac.dH.error)

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
