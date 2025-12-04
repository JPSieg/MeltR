#'Calculates extinction coefficients
#'
#'Calculates the extinction coefficient for a DNA or RNA sequences using the
#'nearest neighbor method. Also calculates the sum of all the extinction
#'coefficients in the input vector. Nearest neighbor parameters are from:
#'
#'Puglisi, J. D.; Tinoco, I. Absorbance Melting Curves of RNA.
#'Methods in Enzymology; RNA Processing Part A: General Methods;
#'Academic Press, 1989; Vol. 180, pp 304-325.
#'https://doi.org/10.1016/0076-6879(89)80108-9.
#'
#'@param NucAcid A vector containing RNA or DNA sequences. Example "c("RNA", "ACCUUUCG", "CGAAAGGU")"
#'@param wavelength The wavelength the data was collected at
#'@param silent To print the extinction coefficient once it is calculated
#'@return A vector containing the extinction coefficient
#' @export
calc.extcoeff = function(NucAcid, wavelength = 260, silent = F){
  tryCatch({
    df.ext = subset(df.ext.data, Wavelength.nm == wavelength) #Grab data for each wavelength
    v_Absorbtivity = df.ext$Absorbtivity.M
    names(v_Absorbtivity) = df.ext$Nucleotide
    l_v_NucAcid = strsplit(NucAcid[-which(NucAcid == NucAcid[1])], split = "") #Makes a list of vectors containing mononucleotides from strings supplied in NucAcid
    Ext = c() #Empty vector of extinction coefficients
    for (i in c(1:length(l_v_NucAcid))){
      #Check strand
      if (silent){}else{
        print("Sequence:")
        print(NucAcid[i+1])
      }
      v_DiNuc = c() #Empty vector of dinucleotides
      v_monoNuc = c() #Empty vector of internal mononucleotides
      v_eDiNuc = c() #Empty vector of dinucleotides
      v_emonoNuc = c() #Empty vector of internal mononucleotides
      for (j in c(1:(length(l_v_NucAcid[[i]])-1))){ #Make a vector of dinucleotides
        DiNuc = paste(l_v_NucAcid[[i]][j], "p", l_v_NucAcid[[i]][(j+1)], sep = "")
        v_eDiNuc = c(v_eDiNuc, DiNuc)
        v_DiNuc = c(v_DiNuc, v_Absorbtivity[[which(names(v_Absorbtivity) == DiNuc)]])
      }
      #Check DiNuc
      if (silent){}else{
        print("Dinucleotides:")
        print(v_eDiNuc)
        print(v_DiNuc)
      }
      for (j in c(2:(length(l_v_NucAcid[[i]])-1))){ #Make a vector of mononucleotides (not including ends)
        MonoNuc = paste(l_v_NucAcid[[i]][j], "p", sep = "")
        v_emonoNuc = c(v_emonoNuc, MonoNuc)
        v_monoNuc = c(v_monoNuc, v_Absorbtivity[[which(names(v_Absorbtivity) == MonoNuc)]])
      }
      #Check monoNuc
      if (silent){}else{
        print("Mononucleotides:")
        print(v_emonoNuc)
        print(v_monoNuc)
      }
      #Subract internal mononucleotides from double the external mononucleotides
      Ext[i] = 2*sum(v_DiNuc) - sum(v_monoNuc)
    }
    #Clean up and format a list of extinction coefficients
    extcoef = list()
    extcoef[[1]] = sum(Ext)
    for (i in c(1:length(Ext))){
      extcoef[[i+1]] = Ext[i]
    }
    names(extcoef) <- c("Total", NucAcid[c(2:length(NucAcid))])
    if (silent){}else{
      print("Final:")
      print(extcoef)
    }
    return(extcoef)
  },
  error = function(e){print("There is no nucleotide extinction coefficient at this wavelength for a nucleotide you specified in your sequence. You will need to provide your own custom extinction coefficient")})
}

