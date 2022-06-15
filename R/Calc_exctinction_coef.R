#'Calculates extinction coefficents
#'
#'Calculates the extinnction coefficient for a DNA or RNA sequences using the
#'nearest neigbor method. Also calculates the sum of all the extinction
#'coefficients in the input vector. Nearest neighbor parameters are from:
#'
#'Puglisi, J. D.; Tinoco, I. Absorbance Melting Curves of RNA.
#'Methods in Enzymology; RNA Processing Part A: General Methods;
#'Academic Press, 1989; Vol. 180, pp 304-325.
#'https://doi.org/10.1016/0076-6879(89)80108-9.
#'
#'@param NucAcid A vector containing RNA or DNA sequences. Example "c("RNA", "ACCUUUCG", "CGAAAGGU")"
#'@return A data frame containing the information in the .ct file
#' @export
calc.extcoeff = function(NucAcid, wavelength = 260){
  tryCatch({

    df.ext = subset(df.ext.data, Wavelength.nm == wavelength)
    RNA = df.ext$Absorbtivity.M
    names(RNA) = df.ext$Nucleotide
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
  },
  error = function(e){print("There is no nucleotide extinction coefficient at this wavelength for a nucleotide you specified in your sequence. You will need to provide your own custom extinction coefficient")})
  output = extcoef
  print(output)
  output <- output
}
