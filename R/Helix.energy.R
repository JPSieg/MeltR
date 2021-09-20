#'Calculates folding INN-HB parameters for a helix sequnece
#'
#'@param seqF Forward RNA sequence
#'@param seqR Reverse RNA sequence
#'@param output The format you want the data to be in. "df" for a dataframe counting the parameters or "energy" for the energy of folding.
#'@param F.Q Is there a fluorophore or a quencher on the helix.
#'@param AA.UU = -0.9, #UU.AA 1 Energy for a folding parameter
#'@param AU.AU = -1.1,#AU.AU Energy for a folding parameter
#'@param UA.UA = -1.3, #UA.UA Energy for a folding parameter
#'@param CU.AG = -2.1, #AG.CU 2 Energy for a folding parameter
#'@param CA.UG = -2.1, #UG.CA 3 Energy for a folding parameter
#'@param GU.AC = -2.2, #AC.GU 4 Energy for a folding parameter
#'@param GA.UC = -2.1, #UC.GA 5 Energy for a folding parameter
#'@param CG.CG = -2.4, #CG.CG Energy for a folding parameter
#'@param GG.CC = -3.3, #CC.GG 6 Energy for a folding parameter
#'@param GC.GC = -3.4, #GC.GC Energy for a folding parameter
#'@param Initiation = 4.09 Energy for a folding parameter
#'@param Term.AU = 0.45 Energy for a folding parameter
#'@param FAMC.GBHQ1 = -3.93 Energy for a folding parameter
#'@param Symmetry = 0.43 Energy for a folding parameter
#'@return A data frame containing the data in the "information in the ".o3a" formatted file.
#' @export
Helix.energy = function(seqF = "UUCCCU",
                       seqR = "AGGGAA",
                       output = "df",
                       F.Q = FALSE,
                       Fluor = "FAM",
                       Qunecher = "BHQ1",
                       AA.UU = -0.9, #UU.AA 1
                       AU.AU = -1.1,#AU.AU
                       UA.UA = -1.3, #UA.UA
                       CU.AG = -2.1, #AG.CU 2
                       CA.UG = -2.1, #UG.CA 3
                       GU.AC = -2.2, #AC.GU 4
                       GA.UC = -2.1, #UC.GA 5
                       CG.CG = -2.4, #CG.CG
                       GG.CC = -3.3, #CC.GG 6
                       GC.GC = -3.4, #GC.GC
                       Initiation = 4.09,
                       Term.AU = 0.45,
                       FAMC.GBHQ1 = -3.93,
                       FAMU.ABHQ1 = -3.0,
                       FAMG.CBHQ1 = -4,
                       FAMA.UBHQ1 = -2.5,
                       Symmetry = 0.43){
  ####Split sequence string####

  vector.seqF = strsplit(seqF, split = "")[[1]]
  vector.seqR = strsplit(seqR, split = "")[[1]]

  ####If Fluorophore/quencher labeled####

  if (F.Q){
    vector.seqF = c("FAM", vector.seqF)
    vector.seqR = c(vector.seqR, "BHQ1")
  }

  ####Fill in NN terms####

  vector.terms = c()

  for (i in 1:(length(vector.seqF) - 1)){
    vector.terms[i] = paste(vector.seqF[i], vector.seqF[i+1], ".", rev(vector.seqR)[i+1], rev(vector.seqR)[i], sep = "")
  }

  ####Add AU term penalty####

  AU.term = FALSE

  if(vector.seqF[length(vector.seqF)] == "A"){
    AU.term = TRUE
  }

  if(vector.seqF[length(vector.seqF)] == "U"){
    AU.term = TRUE
  }

  if (AU.term){
    vector.terms = c("Term.AU", vector.terms)
  }

  AU.start = FALSE

  if(vector.seqF[1] == "A"){
    AU.start = TRUE
  }

  if(vector.seqF[1] == "U"){
    AU.start = TRUE
  }

  if (AU.start){
    vector.terms = c("Term.AU", vector.terms)
  }

  ####Add initiation penalty####

  vector.terms = c("Initiation", vector.terms)

  ####Add symmetry penalty####

  if (seqF == seqR){
    vector.terms = c("Symmetry", vector.terms)
  }

  ####Filter out redundant params####

  for (i in 1:length(vector.terms)){
    if(vector.terms[i] == "UU.AA"){ #1
      vector.terms[i] = "AA.UU"
    }
    if(vector.terms[i] == "AG.CU"){ #2
      vector.terms[i] = "CU.AG"
    }
    if(vector.terms[i] == "UG.CA"){ #3
      vector.terms[i] = "CA.UG"
    }
    if(vector.terms[i] == "AC.GU"){ #4
      vector.terms[i] = "GU.AC"
    }
    if(vector.terms[i] == "UC.GA"){ #5
      vector.terms[i] = "GA.UC"
    }
    if(vector.terms[i] == "CC.GG"){ #5
      vector.terms[i] = "GG.CC"
    }
  }

  ####Generate equation and parse####

  for (i in 2:length(vector.terms)){
    if (i == 2){
      string.formula = paste(vector.terms[i-1], vector.terms[i], sep = " + ")
    }else{
      string.formula = paste(string.formula, vector.terms[i], sep = " + ")
    }
  }
  df.param = data.frame(AA.UU = length(which(vector.terms == "AA.UU")), #UU.AA 1
                        AU.AU = length(which(vector.terms == "AU.AU")), #AU.AU
                        UA.UA = length(which(vector.terms == "UA.UA")), #UA.UA
                        CU.AG = length(which(vector.terms == "CU.AG")), #AG.CU 2
                        CA.UG = length(which(vector.terms == "CA.UG")), #UG.CA 3
                        GU.AC = length(which(vector.terms == "GU.AC")), #AC.GU 4
                        GA.UC = length(which(vector.terms == "GA.UC")), #UC.GA 5
                        CG.CG = length(which(vector.terms == "CG.CG")), #CG.CG
                        GG.CC = length(which(vector.terms == "GG.CC")), #CC.GG 6
                        GC.GC = length(which(vector.terms == "GC.GC")), #GC.GC
                        Term.AU = length(which(vector.terms == "Term.AU")),
                        FAMC.GBHQ1 = length(which(vector.terms == "FAMC.GBHQ1")),
                        FAMU.ABHQ1 = length(which(vector.terms == "FAMU.ABHQ1")),
                        FAMG.CBHQ1 = length(which(vector.terms == "FAMG.CBHQ1")),
                        FAMA.UBHQ1 = length(which(vector.terms == "FAMA.UBHQ1")),
                        Symmetry = length(which(vector.terms == "Symmetry")),
                        Initiation = 1)

  fold.E = eval(parse(text = string.formula))

  print(fold.E)

  if (output == "df"){
    output = df.param
  }else{
    output = fold.E
  }
  output = output
}
