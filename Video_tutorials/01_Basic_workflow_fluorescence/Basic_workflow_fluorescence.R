
####Find the right spot in the computer####

setwd("~/Jacob/R_packages/MeltR/Video_tutorials/01_Basic_workflow_fluorescence")

####Load in data####

list.files()
df = read.csv("js5060_Helix_J_formatted.csv")

####Load packages####

library(tidyverse)
library(ggrepel)
library(MeltR)

####Check data####

head(df)

ggplot(df %>% filter(Reading == 1),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel() +
  theme_classic()

ggplot(df %>% filter(Reading == 30),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel() +
  theme_classic()

####Fit data####

?meltR.F
meltR.F(df,
        Save_results = "all")

####Refine####

meltR.F(df,
        Save_results = "all",
        Kd_range = c(50, 500),
        Kd_error_quantile = 0.4)
