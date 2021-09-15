#'Enumerates all the possible base pairing states for a RNA helix and compiles a list of thermodynamic parameters for each state
#'
#'Generates all possible base pairing states for a helix. Assumes that strands are not allowed to slips.
#'
#'@param NucAcid A vector containing RNA or DNA sequences. Flourophores and quenchers should be separated from the RNA sequnece with a dash. Example "c("RNA", "FAM-ACCUUUCG", "CGAAAGGU-BHQ1")"
#'@return A list of data frames containing the sequence of the A strand, the sequence of the B strand, and a 1 if a BP is present or a 0 if no BP is present.e
#' @export
enumerate.states = function(NucAcid = c("RNA", "F-CGAAAGGU","ACCUUUCG-Q")){
  #NucAcid = c("RNA", "F-ACCUUUCG","CGAAAGGU-Q")

  A = c(strsplit(NucAcid[2], split = "-")[[1]][1], strsplit(strsplit(NucAcid[2], split = "-")[[1]][2], split = "")[[1]])
  B = c(strsplit(NucAcid[3], split = "-")[[1]][2], rev(strsplit(strsplit(NucAcid[3], split = "-")[[1]][1], split = "")[[1]]))

  Term.AU = 0

  INN.HB.P = c()

  for (i in 1:(length(A) - 1)){
    INN.HB.P[i] = paste(A[i], A[i+1], "/", B[i+1], B[i], sep = "")
  }

  df.INN.HB = data.frame("INN.HB.Pairs" =  INN.HB.P)

  matrix.states <- data.matrix(prob::tosscoin(length(df.INN.HB$INN.HB.Pairs)))

  matrix.states[which(matrix.states == 2)] <- 0

  list.df.INN.HB <- {}
  list.formula.state <- {}
  na = c()

  for (i in 1:nrow(matrix.states)){
    list.df.INN.HB[[i]] <- df.INN.HB
    list.df.INN.HB[[i]]$BP <- as.numeric(matrix.states[i,])
    a <- list.df.INN.HB[[i]]
    loop.size = 0
    form = c()
    kill = FALSE
    for (j in 1:(length(a$INN.HB.Pairs))){
      if (kill){
        #list.df.INN.HB[[i]] <- NA
        form <- NA
      }else{
        if (loop.size == 0){
          if (a$BP[j] == 1){
            form = c(form, as.character(a$INN.HB.Pairs[j]))
          }else{
            loop.size = loop.size + 1
            if (j - length(a$INN.HB.Pairs) != 0){
              if (a$BP[j + 1] == 1){
                if (j != 1){
                  kill = TRUE
                }
              }
            }
          }
        }else{
          if (a$BP[j] == 1){
            if (is.null(form)){
              form = c(form, as.character(a$INN.HB.Pairs[j]))
              loop.size = 0
            }else{
              form = c(form, paste(loop.size - 1, "N loop"), as.character(a$INN.HB.Pairs[j]))
              loop.size = 0
              if (loop.size - 1 == 0){kill <- TRUE}
            }
          }else{
            loop.size = loop.size + 1
          }
        }
      }
    }
    if (length(form) != 0){
      if (length(form) > 1){
        if(length(a$BP[which(a$BP == 0)]) != length(a$BP)){form = c(form, "Initiation")}
        if (A[length(A)] == "A"){form = c(form, "Term AU penalty")}
        if (A[length(A)] == "U"){form = c(form, "Term AU penalty")}
        list.formula.state[[i]] <- form
      }else{
        if(is.na(form)){
          na = c(na, i)
          }else{
          if(length(a$BP[which(a$BP == 0)]) != length(a$BP)){form = c(form, "Initiation")}
          if (A[length(A)] == "A"){form = c(form, "Term AU penalty")}
          if (A[length(A)] == "U"){form = c(form, "Term AU penalty")}
            list.formula.state[[i]] <- form
        }
      }
    }else{
      form = "Unpaired"
      list.formula.state[[i]] <- form
    }

  }

  list.df.INN.HB <- list.df.INN.HB[-na]
  list.formula.state <- list.formula.state[-na]

  output <- list(list.df.INN.HB, list.formula.state, NucAcid)
  names(output) <- c("States", "Formulas", "Sequence")
  output <- output

}

