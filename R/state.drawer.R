#'Draws all possible states for a specified double helix
#'
#'Generates R2R inputs for all states in an enumerate.states output object. Run "r2r --disable-usage-warning" on the master R2R meta file in order
#'to generate secondary structural depictions of all states.
#'
#'@param States A MeltR enumerate.states output for a helix.
#' @return Writes a directory of R2R stockholm files and a master R2R meta file you can use to draw all of the states at once.
state.drawer = function(states = enumerate.states()){
  ####Compile list of names for wrtiting files and directories based on the sequence####
  Nuc_acid.name = paste(states$Sequence[1], states$Sequence[2], states$Sequence[3], sep = "_")

  dir.name = paste(states$Sequence[1], states$Sequence[2], states$Sequence[3], "R2R", "files", sep = "_")

  ####Make a folder to contain the R2R inputs####

  dir.create(dir.name)

  ####Generate a R2R formatted stockholm file for each state####

  paths = c()

  for (i in 1:length(states$States)){
    a = states$States[[i]]
    form = states$Formulas[[i]]
    NucAcid = states$Sequence
    A = c(strsplit(NucAcid[2], split = "-")[[1]][1], strsplit(strsplit(NucAcid[2], split = "-")[[1]][2], split = "")[[1]])
    B = c(strsplit(NucAcid[3], split = "-")[[1]][2], rev(strsplit(strsplit(NucAcid[3], split = "-")[[1]][1], split = "")[[1]]))

    Nucleotide = c(A, "U", "U", "U", "U", rev(B))

    Dotbracket = rep(NA, length(Nucleotide))

    for (j in 1:length(a$BP)) {
      if (a$BP[j] == 1){
        Dotbracket[j] <- "<"
        Dotbracket[j +1] <- "<"
        Dotbracket[length(Dotbracket) - j + 1] <- ">"
        Dotbracket[length(Dotbracket) - j] <- ">"
      }
    }

    Dotbracket[is.na(Dotbracket)] <- "."

    df.r2easyR = data.frame(Nucleotide, Dotbracket)

    paths[i] <- paste(dir.name, "/", Nuc_acid.name, "_state_", i, sep = "")

    R2easyR::r2easyR.write(paths[i], df.r2easyR)

    con = file(paste(dir.name, "/", Nuc_acid.name, "_state_", i, ".sto", sep = ""))

    Lines = readLines(con)

    Lines[8] <- "#=GF R2R tick_label_disable_default_numbering"

    writeLines(Lines, con)

    close(con)

    list.files("RNA_F-ACCUUUCG_CGAAAGGU-Q_R2R_files")

    R2easyR::r2easyR.grey_letters_editor(paste(dir.name, "/", Nuc_acid.name, "_state_", i, ".sto", sep = ""),
                                         c((length(A)+1):(length(A)+4)),
                                         "white")
  }

  ####Write a master R2R meta file####

  Lines <- paste(paths, ".sto\toneseq\tdefault", sep = "")

  con = file(paste(Nuc_acid.name, ".r2r_meta", sep = ""))
  writeLines(Lines, con)
  close(con)
}


