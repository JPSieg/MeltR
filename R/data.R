
# This goes in R/data.R

#' @title df.fluor.data
#' @description Fluorecence binding isotherms measuring FAM-CGUAUGUA binding to UACAUACG-BHQ1 in 240 mM NaCl 140 mM KCl 2 mM MgCl2 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0
#' @format A data frame with 4114 rows and 8 variables:
#' \describe{
#'   \item{\code{Well}}{double Well location in a microplate}
#'   \item{\code{Reading}}{double Plate reading on the same sample}
#'   \item{\code{Temperature}}{double Temperature in degrees Celcius}
#'   \item{\code{A}}{double total FAM-CGUAUGUA concentration in nM}
#'   \item{\code{B}}{double total UACAUACG-BHQ1 concentration in nM}
#'   \item{\code{Emission}}{double Fluorecence emission}
#'   \item{\code{Helix}}{double Helix ID}
#'   \item{\code{Conditition}}{double Buffer condition in the description}
#'}
"df.fluor.data"

# This goes in R/data.R

#' @title df.abs.data
#' @description Absorbance melting curves for a heteroduplex helix composed of CGAAAGGU/ACCUUUCG. Samples contain strands at equamolar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0.
#' @format A data frame with 1800 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{double The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{double The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{double The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{double The absorbance at 260 nm}
#'}
"df.abs.data"

