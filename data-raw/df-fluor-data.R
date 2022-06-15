## code to prepare `df.fluor.data` dataset goes here

df.ext.data = read.csv("data/Ext_coef_data.csv")

usethis::use_data(df.ext.data, overwrite = TRUE)
