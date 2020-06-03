#'Pulls and combines data from a single wavelength from multiple ".o3a" files
#'
#'Reads all ".o3a" files in a directory, pulls out the data from a single wavlength,
#'and combines the data into a data frame for plotting and subsequent data analysis.
#'Also includes an option to normalize the data to a reading and write the data frame
#'to a ".csv" file.
#'
#'@param file Path to a directory containing ".o3a" formatted files.
#'@param write Optional prefix for the ".csv" formatted file you want to write. Example is "Name". Default = FALSE to not write a file.
#'@param Wavelength Wavelength in nm that you want to pull out of the ".o3a" files. Default = "260".
#'@param Remove_readings Optional reading to remove from the data set.
#'@param Pathlength Pathlength in cm use to collect the data in the ".o3a" files, formatted into a vector. Example: Pathlength = c(1, 1, 0.5, 0.5, 0.1). Default = 1 cm.
#'@param Sample_names Optional names of the samples in each ".o3a" file, Formatted into a vector. Example: Sample_names = c("Leo", "Theo", "Betty", "Apallo") or  Sample_names = c(1:4). Default = NULL uses the ".o3a" file names as sample names
#'@param Norm_reading Optional reading to normalize the data to. Example: Norm_reading = 1. Default = FALSE to not normalize the data.
#'@return A data frame.
#' @export
o3a.to.MeltR.csv = function(file,
                            write = FALSE,
                            Wavelength = "260",
                            Remove_readings = c(),
                            Pathlength = rep(1, 1000),
                            Sample_names = "NULL",
                            Norm_reading = FALSE){
  if (Sample_names[1] == "NULL"){
    Sample_names = list.files(file)
  }
  data <- {}
  for (i in c(1:length(list.files(file)))){
    print(list.files(file)[i])
    data[[i]] <- t(read.delim(paste(file, list.files(file)[i], sep = "/"), header = FALSE, row.names = 1))
    colnames(data[[i]])[1] <- "Temperature"
    data[[i]] <- data.frame(data[[i]])
    print(head(data[[i]]))
  }
  output <- data.frame("Sample" = c(), "Pathlength" = c(), "Temperature" = c(), "Absorbance" = c())
  for (i in c(1:length(list.files(file)))){
    a <- {}
    a <- data.frame("Sample" = c(Sample_names[i]),
                    "Pathlength" = c(Pathlength[i]),
                    "Temperature" = c(data[[i]]$Temperature)[-Remove_readings],
                    "Absorbance" = c(data[[i]][ , which(colnames(data[[i]]) == paste("X", Wavelength, sep = "")) ])[-Remove_readings]
    )
    if (Norm_reading != FALSE){
      a$Absorbance <- a$Absorbance/a$Absorbance[Norm_reading]
    }
    output <- rbind(output, a)
  }
  if (write != FALSE){
    write.table(output, paste(write, ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print("Output")
  print(head(output))
  output <- output
}
