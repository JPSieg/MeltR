
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
#' @description Absorbance melting curves for a heteroduplex helix composed of CGAAAGGU/ACCUUUCG. Samples contain strands at equal-molar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0.
#' @format A data frame with 1800 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'}
"df.abs.data"

# This goes in R/data.R

#' @title df.abs.non2S
#' @description Absorbance melting curves for a heteroduplex helix composed of CGCGUUAUAU/AUAUUUCGCG. Samples contain strands at equal-molar concentrations in 1 M NaCl 20 mM MOPS 0.01 mM EDTA 0.001% SDS pH 7.0.
#' @format A data frame with 1448 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'}
"df.abs.non2S"

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

# This goes in R/data.R

#' @title df.test.data
#' @description Modeled absorbance melting curves for a data sets used in the testing protocols.
#' @format A data frame with 14661 rows and 5 variables:
#' \describe{
#'   \item{\code{Sample}}{ The sample melted in different cuvettes. Sample 1 is a buffer blank and other samples contain different concentrations of RNA}
#'   \item{\code{Temperature}}{ The temperature where the absorbance was recorded, in degrees Celcius}
#'   \item{\code{Pathlength}}{ The pathlength of the cuvette for each sample in cm}
#'   \item{\code{Absorbance}}{ The absorbance at 260 nm}
#'   \item{\code{DS}}{ The data set. DS1: A monomolecular 2 state helix with no blanks.  DS2: A Heteroduplex 2 state helix with no blanks. DS3: A Homoduplex 2 state helix with no blanks. DS4: A Monomolecular 2 state helix with 1 blank. DS5: A Heteroduplex 2 state helix with 1 blank. DS6: A Homoduplex 2 state helix with 1 blank. DS7: A Monomolecular 2 state helix with 3 blanks. DS8: A Heteroduplex 2 state helix with 3 blanks.  DS9: A Homoduplex 2 state helix with 3 blanks.}
#'}
"df.test.data"

# This goes in R/data.R

#' @title df.test.parameters
#' @description A data frame describing variables to run the meltR.A testing protocol.
#' @format A data frame with 3888 rows and 9 variables:
#' \describe{
#'   \item{\code{blank}}{The method used to blank the data set}
#'   \item{\code{methods}}{The methods to use in the data set}
#'   \item{\code{Mmodel}}{The molecular model}
#'   \item{\code{NucAcid}}{The nucleic acid or concentrations}
#'   \item{\code{outliers}}{The outliers to remove}
#'   \item{\code{Tm_method}}{The method to determine the Tm in method 2}
#'   \item{\code{Wavelength}}{The wavelength to use for extinction coefficient calculation}
#'   \item{\code{fitTs}}{What baseline trimming to use}
#'   \item{\code{DS}}{The data set to use}
#'}
"df.test.parameters"

# This goes in R/data.R

#' @title df.test.results
#' @description The expected results from the testing protocol if everything is working.
#' @format A data frame with 10368 rows and 9 variables:
#' \describe{
#'   \item{\code{Method}}{The meltR.A method}
#'   \item{\code{dH}}{The enthalpy in kcal/mol}
#'   \item{\code{SE.dH}}{The error in the enthalpy in kcal/mol}
#'   \item{\code{dS}}{The entropy in cal/mol/K}
#'   \item{\code{SE.dS}}{The error in the entropy in cal/mol/K}
#'   \item{\code{dG}}{The Gibbs free energy at 37C in kcal/mol}
#'   \item{\code{SE.dG}}{The error in the Gibbs free energy at 37C in kcal/mol}
#'   \item{\code{Tm_at_0.1mM}}{An expectation maximized Tm at Ct = 0.1 mM in Celcius}
#'   \item{\code{SE.Tm_at_0.1mM}}{The error in the expectation maximized Tm at Ct = 0.1 mM in Celcius}
#'}
"df.test.results"

# This goes in R/data.R

#' @title df.BLtrimmer.test.parameters
#' @description A data frame describing variables to run the BLtrimmer testing protocol.
#' @format A data frame with 30 rows and 3 variables:
#' \describe{
#'   \item{\code{Trim.method}}{The method used to generate trimmed baselines}
#'   \item{\code{Assess.method}}{The method used to assess baseline combinations}
#'   \item{\code{quantile.threshold}}{The quantile threshold for baseline combinations to be passed to the results}
#'}
"df.test.parameters"

# This goes in R/data.R

#' @title df.BLtrimmer.test.results
#' @description The expected results from the testing protocol if everything is working.
#' @format A data frame with 90 rows and 5 variables:
#' \describe{
#'   \item{\code{Method}}{The meltR.A method}
#'   \item{\code{dH}}{The enthalpy in kcal/mol}
#'   \item{\code{dS}}{The entropy in cal/mol/K}
#'   \item{\code{dG}}{The Gibbs free energy at 37C in kcal/mol}
#'   \item{\code{Tm_at_0.1mM}}{An expectation maximized Tm at Ct = 0.1 mM in Celcius}
#'}
"df.BLtrimmer.test.results"
