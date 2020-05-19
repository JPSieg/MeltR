#'Reads an ".o3a" file into R
#'
#'Reads an ".o3a" file, produced by the OLIS absorbance spectrophotometer, into a R
#'data frame. The temperature data is in the first column of the data frame.
#'The absorbance for each wavelength is stored in columns that named by the absorbance
#'in nm. Example, "x260".
#'
#'@param file Path to a single ".o3a" formatted file.
#'@return A data frame containing the data in the "information in the ".o3a" formatted file.
#' @export
read.o3a = function(file){
  output <- t(read.delim(file, header = FALSE, row.names = 1))
  colnames(output)[1] <- "Temperature"
  output <- data.frame(output)
  print(output)
}
