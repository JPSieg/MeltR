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
calc.extcoeff = function(NucAcid){
  RNA <- list(15340, 7600, 12160, 10210, 13650, 10670, 12790, 12140, 10670, 7520, 9390, 8370, 12920, 9190, 11430, 10960, 12520, 8900, 10400, 10110)
  names(RNA) <- c("Ap", "Cp", "Gp", "Up", "ApA", "ApC", "ApG", "ApU", "CpA", "CpC", "CpG", "CpU", "GpA", "GpC", "GpG", "GpU", "UpA", "UpC", "UpG", "UpU")
  DNA <- list(15340, 7600, 12160, 8700, 13650, 10670, 12790, 11420, 10670, 7520, 9390, 7660, 12920, 9190, 11430, 10220, 11780, 8150, 9700, 8610)
  names(DNA) <- c("Ap", "Cp", "Gp", "Tp", "ApA", "ApC", "ApG", "ApT", "CpA", "CpC", "CpG", "CpT", "GpA", "GpC", "GpG", "GpT", "TpA", "TpC", "TpG", "TpT")
  if (NucAcid[1] == "RNA"){
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
    output <- extcoef
  }
  if (NucAcid[1] == "DNA"){
    b <- c()
    b <- strsplit(NucAcid[-which(NucAcid == NucAcid[1])], split = "")
    c <- {}
    d <- {}
    e <- c()
    for (i in c(1:length(b))){
      c[[i]] <- c(1:(length(b[[i]])-1))
      d[[i]] <- c(1:(length(b[[i]])))
      for (j in c(1:(length(b[[i]])-1))){
        c[[i]][j] <- RNA[[which(names(DNA) == paste(b[[i]][j], "p", b[[i]][j+1], sep = ""))]]
      }
      for (j in c(1:(length(b[[i]])))){
        d[[i]][j] <- RNA[[which(names(DNA) == paste(b[[i]][j], "p", sep = ""))]]
      }
      e[i] <- 2*sum(c[[i]]) - sum(d[[i]])
    }
    extcoef <- list()
    extcoef[[1]] <- sum(e)
    for (i in c(1:length(b))){
      extcoef[[i+1]] <- e[i]
    }
    names(extcoef) <- c("Total", NucAcid[c(2:length(NucAcid))])
    output <- extcoef
  }
  print(output)
  output <- output
}
