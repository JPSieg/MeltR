## code to prepare `df.test.parameters` dataset goes here

df.test.parameters = read.csv("data/Parameters_for_test.csv")

usethis::use_data(df.test.parameters, overwrite = TRUE)

## code to prepare `df.BLtrimmer.test.parameters` dataset goes here

df.BLtrimmer.test.parameters = read.csv("data/BLtrimmer_parameters_for_test.csv")

usethis::use_data(df.BLtrimmer.test.parameters, overwrite = TRUE)
