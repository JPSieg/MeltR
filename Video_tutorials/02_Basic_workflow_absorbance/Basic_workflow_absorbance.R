####Find the right place in the computer####

setwd("~/Jacob/R_packages/MeltR/Video_tutorials/02_Basic_workflow_absorbance")

####Load packages####

library(tidyverse)
library(MeltR)

####Read in data####

list.files()

df = read.csv("Heteroduplex_js3041_HelixA_No_label.csv")

####Check data####

head(df)
ggplot(df,
       aes(x = Temperature, y = Absorbance, color = factor(Sample))) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = c(25, 80))

####Fit data####

?meltR.A

fit = meltR.A(df,
              blank = 1,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              fitTs = c(25, 80),
              Save_results = "all")

####Auto-trimming data####

BLTrimmer(fit, Save_results = "all")
