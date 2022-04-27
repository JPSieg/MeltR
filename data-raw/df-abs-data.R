## code to prepare `df.abs.data` dataset goes here

df.abs.data = read.csv("data/Heteroduplex_js3041_HelixA_No_label.csv")

usethis::use_data(df.abs.data, overwrite = TRUE)
