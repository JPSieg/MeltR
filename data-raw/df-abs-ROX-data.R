## code to prepare `df.abs.ROX.data` dataset goes here

df.abs.ROX.data = read.csv("data/js3047_ROX_data.csv")

usethis::use_data(df.abs.ROX.data, overwrite = TRUE)
