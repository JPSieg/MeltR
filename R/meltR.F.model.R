#'Fit fluorescence binding isotherms to obtain thermodynamic parameters
#'
#'Models data for a set of fluorescence binding isotherms given a enthalpy and entropy of folding. Assumes three replicates
#'of quencher labeled strand concentration =  0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, and 1000 nM with
#'a constant concentration of 200 nM FAM. Generates isotherms from 20 to 80 degrees celcius with 0.5 degrees celcius
#'steps. Models error using the rnorm function. 
#'
#'@param H The enthalpy of folding in kcal/mol.
#'@param S The entropy of folding in kcal/mol/K.
#'@param fmax The approximate fluorescence maximum of the fluorophore labeled strand in its single stranded state
#'@param fmin The approximate fluorescence minimum for the double stranded state where the fluorophore labeled strand is base paired with the quencher labeled strand.
#'@param Emission_SD The approximate experimental error for each emission reading.
#'@return A dataframe containing the modeled data. Ready to feed into meltR.F.
#' @export
meltR.F.model = function(H,
                         S,
                         fmax = 2,
                         fmin = 0.2,
                         Emission_SD = 0.05){
  df.model <- data.frame("H" = H,
                         "S" = S,
                         "fmax" = fmax,
                         "fmin" = fmin,
                         "A" = 200,
                         "B" = rep(c(c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                                     c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                                     c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000)), length(seq(20, 80, by = 0.5))),
                         "Well" = rep(c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
                                        "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
                                        "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"), length(seq(20, 80, by = 0.5))))
  
  Temperature <- c()
  Reading <- c()
  
  for (i in 1:length(seq(20, 80, by = 0.5))){
    Temperature <- c(Temperature,
                     rep(seq(20, 80, by = 0.5)[i],
                         length(c(c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                                  c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                                  c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000)))))
    Reading <- c(Reading,
                 rep(c(1:length(seq(20, 80, by = 0.5)))[i],
                     length(c(c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                              c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000),
                              c(0, 1, 10, 50, 100, 150, 200, 250, 400, 600, 800, 1000)))))
  }
  
  df.model$Temperature <- Temperature
  df.model$Reading <- Reading
  
  Emission <- Global(df.model$H, df.model$S, df.model$fmax, df.model$fmin,
                     df.model$Reading, df.model$A, df.model$B, Temperature = df.model$Temperature)
  
  Emission <- Emission + rnorm(length(Emission), 0, Emission_SD)
  
  df.model$Emission <- Emission
  
  output <- df.model
}
