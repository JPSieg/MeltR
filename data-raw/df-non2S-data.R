## code to prepare `df.abs.non2S` dataset goes here

df.abs.non2S = read.csv("data/js6040_multistate_data.csv")

usethis::use_data(df.abs.non2S, overwrite = TRUE)
