
df.set = df.test.data
df = df.test.parameters

####Fit on a loop####

list.fits = {}

pb = txtProgressBar(min = 1, max = nrow(df), initial = 1, style = 3)

for (i in 1:nrow(df)){

  setTxtProgressBar(pb, i)

  #print(paste(i, nrow(df), sep = "/"))
  df.DS = subset(df.set, DS == df$DS[i])

  #Set blank

  if (df$blank[i] == "none"){blank = "none"}
  if (df$blank[i] == "1"){blank = 1}
  if (df$blank[i] == "list"){blank = list(c(2, 1),
                                          c(3, 1),
                                          c(5, 4),
                                          c(6, 4),
                                          c(8, 7),
                                          c(9, 7))}
  #Set nucleic acids

  if (df$NucAcid[i] == "Concentration"){
    if (df$DS[i] == "DS1"){c("Concentration", 3.98e-06, 7.38e-06, 9.03e-06, 1.82e-05, 3.17e-05, 5.70e-05, 1.03e-04, 1.50e-04, 2.34e-04)}
    if (df$DS[i] == "DS2"){c("Concentration", 4.23e-06, 5.50e-06, 1.15e-05, 2.01e-05, 2.79e-05, 2.70e-05, 8.88e-05, 1.57e-04, 2.27e-04)}
    if (df$DS[i] == "DS3"){c("Concentration", 4.64e-06, 6.95e-06, 1.28e-05, 1.85e-05, 2.96e-05, 4.80e-05, 9.24e-05, 1.46e-04, 2.17e-04)}
    if (df$DS[i] == "DS4"){c("Concentration", 7.38e-06, 9.03e-06, 1.82e-05, 3.17e-05, 5.70e-05, 1.03e-04, 1.50e-04, 2.34e-04)}
    if (df$DS[i] == "DS5"){c("Concentration", 7.38e-06, 9.03e-06, 1.82e-05, 3.17e-05, 5.70e-05, 1.03e-04, 1.50e-04, 2.34e-04)}
    if (df$DS[i] == "DS6"){c("Concentration", 6.95e-06, 1.28e-05, 1.85e-05, 2.96e-05, 4.80e-05, 9.24e-05, 1.46e-04, 2.17e-04)}
    if (df$DS[i] == "DS7"){c("Concentration", 7.38e-06, 9.03e-06, 3.17e-05, 5.70e-05, 1.50e-04, 2.34e-04)}
    if (df$DS[i] == "DS8"){c("Concentration", 5.50e-06, 1.15e-05, 2.79e-05, 2.70e-05, 1.57e-04, 2.27e-04)}
    if (df$DS[i] == "DS9"){c("Concentration", 6.95e-06, 1.28e-05, 2.96e-05, 4.80e-05, 1.46e-04, 2.17e-04)}
  }else{
    A = eval(parse(text = df$NucAcid[i]))
  }

  #Set T range

  if (df$fitTs[i] == "none"){
    B = NULL
  }
  if (df$fitTs[i] == "fixed"){
    B = c(15, 95)
  }
  if (df$fitTs[i] == "list"){
    B = list(c(15, 95), #1
             c(16, 94), #2
             c(18, 92), #3
             c(10, 94), #4
             c(6, 95), #5
             c(4, 96), #6
             c(8, 93), #7
             c(14, 94), #8
             c(20, 95)) #9
    if (df$blank[i] == "1"){B = B[2:9]}
    if (df$blank[i] == "list"){B = B[4:9]}
  }

  fit = meltR.A(data_frame = df.DS,
                blank = blank,
                Mmodel = df$Mmodel[i],
                NucAcid = A,
                outliers = df$outliers[i],
                Tm_method = df$Tm_method[i],
                wavelength = df$Wavelength[i],
                fitTs = B,
                Silent = T)

  list.fits[[i]] = fit$Summary

}

close(pb)

####Write results####

df = list.fits[[1]]
for (i in 2:length(list.fits)){
  df = rbind(df, list.fits[[i]])
}

####Test if you got the right result####

testthat::test_that("meltR.A works",
          {
            testthat::expect_equal(df, df.test.results)
          })
