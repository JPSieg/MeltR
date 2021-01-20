#'Takes the first and second derivatives of a melting curves
#'
#'Takes the first and second derivatices of all the melting curves in a meltR data frame. Returns a list of dataframes, one containing the first deeivative and one containing the second.
#'
#'@param data_frame The MeltR data frame you want to take the derivative of. Columns must be named "Sample", "Pathlength", "Temperature", and "Absorbance.
#'@param n The number of data point you want to smooth the first derivative to. Default n = 4.
#'@return A list of data frames. One data frame contains the first derivative and one data frame contains the second derivative.
#' @export
meltR.A.derivative = function(data_frame, n = 4){
  ####First derivative
  a <- {}
  b <- {}
  for (i in 1:length(unique(data_frame$Sample))){
    a <- subset(data_frame, Sample == unique(data_frame$Sample)[i])
    Sample <- unique(data_frame$Sample)[i]
    x <- c()
    y <- c()
    for (j in n:length(a$Sample)){
      y[j] <- (a$Absorbance[j] - a$Absorbance[(j-(n-1))])/(a$Temperature[j] - a$Temperature[j-(n-1)])
      x[j] <- mean(a$Temperature[(j - n + 1):j])
    }
    b[[i]] <- data.frame("Sample" = Sample,
                                    "Temperature" = x,
                                    "First.derivative" = y)
  }
  first <- b[[1]]
  for (i in 2:length(b)){
    first <- rbind(first, b[[i]])
  }
  first <- first[-which(is.na(first$Temperature)),]
  ####Second derivative###
  a <- {}
  b <- {}
  for (i in 1:length(unique(first$Sample))){
    a <- subset(first, Sample == unique(first$Sample)[i])
    Sample <- unique(first$Sample)[i]
    x <- c()
    y <- c()
    for (j in n:length(a$Sample)){
      y[j] <- (a$First.derivative[j] - a$First.derivative[(j-(n-1))])/(a$Temperature[j] - a$Temperature[j-(n-1)])
      x[j] <- mean(a$Temperature[j:(n-1)])
    }
    b[[i]] <- data.frame("Sample" = Sample,
                         "Temperature" = x,
                         "Second.derivative" = y)
  }
  second <- b[[1]]
  for (i in 2:length(b)){
    second <- rbind(second, b[[i]])
  }
  second <- second[-which(is.na(second$Temperature)),]
  ####Assemble results####
  output <- list(first, second)
  names(output) <- c("first", "second")
  print(head(output$first))
  print(head(output$second))
  output <- output
}
