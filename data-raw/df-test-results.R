## code to prepare `df.test.results` data set goes here

df.test.results = read.csv("data/Test_results.csv")

usethis::use_data(df.test.results, overwrite = TRUE)

## code to prepare `df.BLtrimmer.test.results` data set goes here

df.BLtrimmer.test.results = read.csv("data/BLtrimmer_test_results.csv")

usethis::use_data(df.BLtrimmer.test.results, overwrite = TRUE)

