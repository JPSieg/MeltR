#'Reformats an Applied Biosystems qPCR export file for fitting with MeltR
#'
#'This function reads an Applied Biosystems (AB) qPCR export file (Text tab delimited), pulls out Temperature and normalized emission data,
#'and formats it into a R data frame that MeltR can fit. Requires a .csv file index describing what reagents
#'are in each well. The index should have columns labeled "well", "A" for the first reagent,
#'and "B" for the second reagent. Also an option to save the result as a .csv file (see output_name)
#'and to apply a Temperature correction function.
#'
#'@param data_file Path to the tab delimited AB raw data export file.
#'@param index_file Path to the .csv formatted index file.
#'@param output_name Optional name of the .csv file you want to write. Default = NA to not write a file.
#'@param Temp_correction Optional function to apply to correct the temperature readings, determined by callibrating the temperature in the qPCR machine. Default = NA to not apply a correction.
#'@return A data frame containing the Normalized emission and temperature data for the wells you specified with the index file.
#' @export
AB.qPCR.export.to.MeltR.csv = function(data_file, index_file, output_name = NA, Temp_correction = NA){
  ####Read in temperature data####
  Emission_wells <- c()
  Emissions <- list()
  Temperature_wells <- c()
  Temperatures <- list()
  con = file(data_file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    #####Pull out temperature data into a list indexed by well names####
    if (line == "[Melt Region Temperature Data]"){
      while(line != ""){
        line = readLines(con, n = 1)
        if(line == ""){
          break
        }
        str = strsplit(toString(line), split = "\t")
        if (str[[1]][1] == "Well"){ #Make a vector containing the readings
          Temp_Reading <- c()
          a <- c(5:length(str[[1]]))
          for (i in c(1:length(a))){
            Temp_Reading[i] <- strsplit(toString(str[[1]][a[i]]), split = " ")[[1]][2]
          }
        }
        if (str[[1]][1] != "Well"){ #Make a list of vectors containing the temperaure in each well
          Temperature_wells <- c(Temperature_wells, str[[1]][2])
          Temperatures <- c(Temperatures, list(as.numeric(str[[1]][5:length(str[[1]])])))
          names(Temperatures) <- Temperature_wells
        }
      }
    }
    ####Pull out Normalized emission data inot a list in indexed by well names####
    if (line == "[Melt Region Normalized Data]"){
      while(line != ""){
        line = readLines(con, n = 1)
        if(line == ""){
          break
        }
        str = strsplit(toString(line), split = "\t")
        if (str[[1]][1] == "Well"){ #Make a vector containing the readings
          Em_Reading <- c()
          a <- c(5:length(str[[1]]))
          for (i in c(1:length(a))){
            Em_Reading[i] <- strsplit(toString(str[[1]][a[i]]), split = " ")[[1]][2]
          }
        }
        if (str[[1]][1] != "Well"){ #Make a list of vectors containing the temperaure in each well
          Emission_wells <- c(Emission_wells, str[[1]][2])
          Emissions <- c(Emissions, list(as.numeric(str[[1]][5:length(str[[1]])])))
          names(Emissions) <- Emission_wells
        }
      }
    }
  }
  close(con)
  index <- read.csv(index_file)
  output <- data.frame("Well" = c(), "Reading" = c(), "Temperature" = c(), "A" = c(), "B" = c(), "Emission" = c())
  for (i in 1:length(index$Well)){
    output <- rbind(output,
                    data.frame("Well" = names(Emissions)[which(names(Emissions) == index$Well[i])],
                               "Reading" = Em_Reading,
                               "Temperature" = Temperatures[[which(names(Temperatures) == index$Well[i])]],
                               "A" = index$A[i],
                               "B" = index$B[i],
                               "Emission" = Emissions[[which(names(Emissions) == index$Well[i])]]))
  }
  if (is.na(Temp_correction) == FALSE){
    output$Temperature <- Temp_correction(output$Temperature)
  }
  print(head(output))
  if (is.na(output_name) == FALSE){
    write.csv(output, output_name, row.names = FALSE)
  }
  output <- output
}
