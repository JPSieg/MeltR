## code to prepare `df.abs.ROX.data` dataset goes here

df.FAM.C.BHQ1.data = read.csv("data/js3042_HelixC_labeled.csv")

usethis::use_data(df.FAM.C.BHQ1.data, overwrite = TRUE)
