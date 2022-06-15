## code to prepare `df.abs.ROX.data` dataset goes here

df.abs.ROX.data = read.csv("data/Ext_coef_data.csv")

usethis::use_data(df.abs.ROX.data, overwrite = TRUE)
