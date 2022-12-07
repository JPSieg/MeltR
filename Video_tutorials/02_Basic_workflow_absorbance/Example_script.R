####Find the right place in the computer####

setwd("~/Jacob/R_packages/MeltR/Video_tutorials/02_Basic_workflow_absorbance")

####load in packages####

library(MeltR)
library(tidyverse)

####Load in the data####

list.files()
df = read.csv("Heteroduplex_js3041_HelixA_No_label.csv")

####Check data####

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance, color = factor(Sample))) +
  geom_point() +
  geom_vline(xintercept = c(25, 80))

####Fit the data with MeltR####

#meltR.A

?meltR.A

meltR.A.fit = meltR.A(df,
        blank = 1,
        NucAcid = c("RNA", "CGAAGGU", "ACCUUUCG"),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(25, 80),
        Save_results = "all")

#BLTrimmer 

?BLTrimmer

Trimmed.fit = BLTrimmer(meltR.A.fit, Save_results = "all")
