
# This goes in R/data.R

#' @title df.fluor.data
#' @description Fluorecence binding isotherms measuring FAM-CGUAUGUA binding to UACAUACG-BHQ1 in 240 mM NaCl 140 mM KCl 2 mM MgCl2 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0
#' @format A data frame with 4114 rows and 8 variables:
#' \describe{
#'   \item{\code{Well}}{ Well location in a microplate}
#'   \item{\code{Reading}}{ Plate reading on the same sample}
#'   \item{\code{Temperature}}{ Temperature in degrees Celcius}
#'   \item{\code{A}}{ total FAM-CGUAUGUA concentration in nM}
#'   \item{\code{B}}{ total UACAUACG-BHQ1 concentration in nM}
#'   \item{\code{Emission}}{ Fluorecence emission}
#'   \item{\code{Helix}}{ Helix ID}
#'   \item{\code{Conditition}}{ Buffer condition in the description}
#'}
"df.fluor.data"

# This goes in R/data.R

#' @title df.abs.data
#' @description Absorbance melting curves for a heteroduplex helix composed of CGAAAGGU/ACCUUUCG. Samples contain strands at equamolar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0.
#' @format A data frame with 1800 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'}
"df.abs.data"

# This goes in R/data.R

#' @title df.abs.ROX.data
#' @description Absorbance melting curves for a heteroduplex helix composed of CGAAAGGU/ACCUUUCG with variable amounts of the ROX fluorophore. Samples contain strands at equamolar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0.
#' @format A data frame with 1800 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'}
"df.abs.ROX.data"

# This goes in R/data.R

#' @title df.ext.data
#' @description Data for calculating extinction coefficients from DNA and RNA sequences
#' @format A data frame with 323 rows and 5 variables:
#' \describe{
#'   \item{\code{Wavelength.nm}}{The wavelength that the extinction coefficient will be calculated}
#'   \item{\code{Absorbtivity.mM}}{The absorbtivity in mM for a nucleotide}
#'   \item{\code{Absorbtivity.M}}{The absorbtivity in M for a nucleotide}
#'   \item{\code{Nucleotide}}{The nucleotide identity}
#'   \item{\code{Reference}}{The reference the absorbtivity values come from}
#'}
"df.ext.data"

# This goes in R/data.R

#' @title df.FAM.C.BHQ1.data
#' @description Absorbance melting curves for a heteroduplex helix composed of FAM-CUGAGUC/GACUCAG-BHQ1. Samples contain strands at equamolar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0. The data has been background subtracted.
#' @format A data frame with 1246 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'}
"df.FAM.C.BHQ1.data"
