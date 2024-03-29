% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calc_exctinction_coef.R
\name{calc.extcoeff}
\alias{calc.extcoeff}
\title{Calculates extinction coefficients}
\usage{
calc.extcoeff(NucAcid, wavelength = 260, silent = F)
}
\arguments{
\item{NucAcid}{A vector containing RNA or DNA sequences. Example "c("RNA", "ACCUUUCG", "CGAAAGGU")"}

\item{wavelength}{The wavelength the data was collected at}

\item{silent}{To print the extinction coefficient once it is calculated}
}
\value{
A vector containing the extinction coefficient
}
\description{
Calculates the extinction coefficient for a DNA or RNA sequences using the
nearest neighbor method. Also calculates the sum of all the extinction
coefficients in the input vector. Nearest neighbor parameters are from:
}
\details{
Puglisi, J. D.; Tinoco, I. Absorbance Melting Curves of RNA.
Methods in Enzymology; RNA Processing Part A: General Methods;
Academic Press, 1989; Vol. 180, pp 304-325.
https://doi.org/10.1016/0076-6879(89)80108-9.
}
