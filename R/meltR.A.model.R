#'Models realistic absorbance melting curves
#'
#'@param NucAcid A vector containing the Nucleic acid type and the sequences you are fitting for calculating extinction coefficients. Examples: c("RNA", "UUUUUU", "AAAAAA"), c("DNA", "GCTAGC"), etc... . For a custom extinction coefficient enter "Custom" followed by the molar extinction coefficients for every nucleic acid in the sample. For example, c("Custom", 10000, 20000).
#'@param Samples The number of samples you want to model
#'@param Absorbance_error The amount of absorbance error
#'@param Mmodel The molecular model you want to fit. Options: "Monomolecular.2State", "Monomolecular.3State", "Heteroduplex.2State", "Homoduplex.2State".
#'@param Tmodel The thermodynamic model you want to fit. Options: "VantHoff". Default = "VantHoff".
#'@return A data frame containing modeled data
#' @export
meltR.A.model = function(NucAcid = c("RNA", "CGCGCG"),
                   Samples = 9,
                   Absorbance_error = 0.01,
                   Mmodel = "Homoduplex.2State",
                   Tmodel = "VantHoff") {
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
      calcTM = function(H, S, Ct){(H/(S/1000)) - 273.15}
      calcTM.SE = function(H, S, SE.H, SE.S, Ct){
        ((H/(S/1000)) - 273.15)*sqrt((SE.H/H)^2 + (SE.S/S)^2)
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
      calcTM = function(H, S, Ct){NA}
      calcTM.SE = function(H, S, SE.H, SE.S, Ct){
        NA
      }
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
      calcTM = function(H, S, Ct){(H/((S/1000) - 0.0019872*log(4/Ct))) - 273.15}
      calcTM.SE = function(H, S, SE.H, SE.S, Ct){
        ((H/((S/1000) - 0.0019872*log(4/Ct))) - 273.15)*sqrt((SE.H/H)^2 + (SE.S/S)^2)
      }
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
      calcTM = function(H, S, Ct){(H/((S/1000) - 0.0019872*log(1/Ct))) - 273.15}
      calcTM.SE = function(H, S, SE.H, SE.S, Ct){
        ((H/((S/1000) - 0.0019872*log(1/Ct))) - 273.15)*sqrt((SE.H/H)^2 + (SE.S/S)^2)
      }
    }
  }

  ####Determine thermodynamic parameters####

  if (Mmodel == "Homoduplex.2State"){
    seqF = NucAcid[2]
    seqR = NucAcid[2]
  }
  if (Mmodel == "Heteroduplex.2State"){
    seqF = NucAcid[2]
    seqR = NucAcid[3]
  }

  print("Predicted dG in kcal/mol")
  dG = MeltR::Helix.energy(seqF, seqR, output = "energy")
  print("Predicted dH in kcal/mol")
  dH = MeltR::Helix.energy(seqF, seqR, output = "energy",
                           AA.UU = -6.82,
                           AU.AU = -9.38,
                           UA.UA = -7.69,
                           CU.AG = -10.48,
                           CA.UG = -10.44,
                           GU.AC = -11.40,
                           GA.UC = -12.44,
                           CG.CG = -10.64,
                           GG.CC = -13.39,
                           GC.GC = -14.88,
                           Initiation = 3.61,
                           Term.AU = 3.72,
                           Symmetry = 0)

  dS = (dH - dG)/(273.15 + 37)

  print("Predicted dS in kcal/mol/K")

  print(dS)

  print("Predicted Tm at 0.1 mM")

  Tm.at.0.1mM = (1/TmModel(dH, dS, log(10^-4))) -273.15
  print(Tm.at.0.1mM)

  ####Generate extinction coefficients####

  Extinct = MeltR::calc.extcoeff(NucAcid)
  Extinct = Extinct$Total

  Extinct.lower = Extinct*0.6


  ####Generate sample concentrations####

  Sample = 1:Samples


  start.Ct = 0.2/Extinct.lower
  end.Ct = 2/(Extinct*0.1)

  lnCt = seq(log(start.Ct), log(end.Ct), length.out = Samples)

  Ct = exp(lnCt)

  print("Concentrations (M)")
  print(Ct)

  ####Generate sample pathlengths####

  A.1cm = Ct*Extinct
  A.0.5cm = Ct*Extinct*0.5

  Pathlength = rep(NA, length(Ct))

  Pathlength[which(A.1cm <= 2)] = 1
  Pathlength[which(A.1cm >= 2)] = 0.5
  Pathlength[which(A.0.5cm >= 2)] = 0.1

  ####Generate slope and intercept####

  mED = c()
  mSS = c()
  bED = c()
  bSS = c()

  for (i in 1:Samples){
    mED[i] = rnorm(1, 0.001113342, 0.0007760626)
    mSS[i] = rnorm(1, 0.001113342, 0.0007760626)
    bED[i] = Ct[i]*Pathlength[i]*Extinct.lower - mED[i]*90
    bSS[i] = Ct[i]*Pathlength[i]*Extinct - mED[i]*90
  }

  ####Consoilidate sample data####

  df = data.frame(seqF, seqR, dG, dH, dS, Tm.at.0.1mM, Extinct, Extinct.lower, mED, mSS, bED, bSS, Sample, Ct,  Pathlength)

  output = list(df)

  ####Model Absorbance data####

  list.df = {}

  for (i in 1:nrow(df)){
    Sample = i
    Temperature = seq(5, 95, 0.5)
    Pathlength = df$Pathlength[i]
    Tm = calcTM(df$dH[i], 1000*df$dS[i], df$Ct[i])
    Absorbance = Model(df$dH[i], Tm, df$mED[i], df$bED[i], df$mSS[i], df$bSS[i], Temperature, df$Ct[i]) + rnorm(length(Temperature), 0, Absorbance_error)
    list.df[[i]] = data.frame(Sample, Temperature, Pathlength, Absorbance)
  }

  df = list.df[[1]]
  i = 2
  while(i < (1 + Sample)){
    df = rbind(df, list.df[[i]])
    i = i + 1
  }

  output[[2]] = df

  output <- output
}
