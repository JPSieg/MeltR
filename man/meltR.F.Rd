% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meltR.f.R
\name{meltR.F}
\alias{meltR.F}
\title{Fit fluorescence binding isotherms to obtain thermodynamic parameters}
\usage{
meltR.F(
  df,
  Kd_error_quantile = 0.25,
  Kd_range = c(10, 1000),
  Start_K = 0.1,
  vh_start = list(H = -70, S = -0.18),
  Optimize_conc = TRUE,
  Low_reading = "auto",
  low_K = FALSE,
  B.conc.Tm = 250,
  Save_results = "none",
  file_prefix = "Fit",
  file_path = getwd(),
  silent = FALSE
)
}
\arguments{
\item{df}{data frame containing fluorescence binding data}

\item{Kd_error_quantile}{Quantile for uncertainty in equilibrium constants for fitting to thermodynamic models. K_error = K standard error/K. Default = 0.25 or the 25\% most accurate binding constants in the Kd_range (See "Kd_range").}

\item{Kd_range}{A custom acceptable nM Kd range to fit for your experiment. Options FALSE for no custom Kd_range or c(start, end) to set a Kd_range. Example: "Kd_range = c(5, 100)"}

\item{Start_K}{A Kd value to start non-linear regression. Default = 0.1.}

\item{vh_start}{A list of initial guesses for the helix folding enthalpy and entropy. By default vh_start = list(H = -70, S = -0.180)}

\item{Optimize_conc}{Deals with a fundamental experimental uncertainty in determination of A =fluorophore and B = Quencher concentrations in the experiment. If TRUE, meltR.f will optimize the concentration for the fluorophore labeled strand based on the shape of low temperature isotherms.}

\item{Low_reading}{Used by the concentration optimization algorithm. The isotherm, or reading, that you want to use to optimize the concentration. Default = "auto" will use the lowest temperature reading.}

\item{low_K}{Used by the concentration optimization algorithm. A low Kd value in nanomolar, that is used to find an optimum ration between A and B strands in the experiment. Default = FALSE to allow the low_K to float in the concentration optimization algorithm. It is best to use values between 0.1 and 10.}

\item{B.conc.Tm}{Only use quencher (or B strands) higher than this threshold in the 1/Tm versus lnCt fitting method, method 3}

\item{Save_results}{What results to save. Options: "all" to save PDF plots and ".csv" formatted tables of parameters, "some" to save ".csv" formatted tables of parameters, or "none" to save nothing.}

\item{file_prefix}{Prefix that you want on the saved files.}

\item{file_path}{Path to the directory you want to save results in.}

\item{silent}{Set to TRUE to run in silent mode (which does not print results in the console). Good for running in loops. Default is TRUE.}
}
\value{
A meltR.F fit object containing a list of data objects containing raw data, data, transformation, fit objects, and statistics from the fits for plotting, exporting, and advanced analysis.
\itemize{
\item 1. VantHoff - A data frame containing the duplex formation energies.
\item 2. K - A data frame containing the results from fitting each isotherm individually.
\item 3. VH_method_1_fit - A nls object containing the fit obtained from the fit obtained from the Van't Hoff plot.
\item 4. VH_method_2_fit - A nls object containing the global fit.
\item 5. Raw_data - The raw data passed back out of MeltR.F with no modifications.
\item 6. First_derivative - The first derivative of each sample. Useful for qualitative comparison of data between conditions..
\item 7. Tms - The approximate Tm of each sample obtained from the maximum of the first derivative. Useful for qualitative comparison of data between solution conditions.
\item 8. R - The mole ratio of fluorophore and quencher labeled RNA, used in the concentration optimization algorithm. The mole ratio “R” was labeled “X” in the theory section to avoid confusion with the gas constant.
\item 9. Fractional_error_between_methods - The amount thermodynamic parameters vary between methods.
}
}
\description{
Automates the trivial but time-consuming tasks associated with non-linear regression.
Uses two non-linear regression methods to calculate thermodynamic parameters. Method 1
fits each isotherm individually then calculates thermodynamic parameters with a Van't Hoff
plot. Method 2 calculates thermodynamic parameters with a global fit, where H and S are constant
between isotherms and the Fmax and Fmin are allowed to float. Also includes an algorithm that
optimizes the mole ratio of fluorophore labeled strands to quencher labeled strands.
}
