####Find the right place in the computer####

setwd("~/Jacob/R_packages/MeltR/Video_tutorials/01_Basic_workflow_fluorescence")

####Load in packages####

library(MeltR)
library(tidyverse)
library(ggrepel)

####Load in the data####

list.files()

df = read.csv("js5060_Helix_J_formatted.csv")

####Look at data####

head(df)

ggplot(df %>% filter(Reading == 100), aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

####Fit data with meltR####

?meltR.F

fit = meltR.F(df, Save_results = "all",
        Kd_error_quantile = 0.48,
        Kd_range = c(50, 500))

fit$Tms
coef(summary(fit$VH_method_1_fit))
