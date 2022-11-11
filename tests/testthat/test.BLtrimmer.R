
####preliminary fit data in MeltR####

fit = meltR.A(df.abs.data,
              blank = 1,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              Silent = T)

df = df.BLtrimmer.test.parameters


####Fit data using all combinations of the important BLTrimmer arguments####

list.df = {}
pb = txtProgressBar(min = 1, max = nrow(df), initial = 1, style = 3)
for (i in 1:nrow(df)){
  setTxtProgressBar(pb, i)
  #print(paste("Test ", i, 30, sep = "/"))
  trim = BLTrimmer(fit,
                   n.combinations = 100,
                   Trim.method = df$Trim.method[i],
                   Assess.method = df$Assess.method[i],
                   quantile.threshold = df$quantile.threshold[i],
                   Silent = T)
  list.df[[i]] = trim$Ensemble.energies
}

####Consolidate into a single data frame####

df = list.df[[1]][,c(1,2,4,6,8)]
for (i in 2:length(list.df)){
  df = rbind(df, list.df[[i]][,c(1,2,4,6,8)])
}

####Determine if results are different####

df.test = df.BLtrimmer.test.results
df.test
df
M.test = as.matrix(df.test[,c(2:5)])
M = as.matrix(df[,c(2:5)])
ME = 100*abs((M - M.test)/M.test)

mean.EdH = mean(ME[,1])
mean.EdS = mean(ME[,2])
mean.EdG = mean(ME[,3])
mean.ETm = mean(ME[,4])

result = c(F, F, F, F)

if (mean.EdH <= 3){result[1] = T}
if (mean.EdS <= 3){result[2] = T}
if (mean.EdG <= 1){result[3] = T}
if (mean.ETm <= 1){result[4] = T}

####Test if results are within the expected error####

testthat::test_that("BLTrimmer works",
                    {
                      testthat::expect_equal(result, c(T,T,T,T))
                    })

