## code to prepare `df.fluor.data` dataset goes here

df.fluor.data = read.csv("data/js5060_Helix_J_formatted.csv")

usethis::use_data(df.fluor.data, overwrite = TRUE)
