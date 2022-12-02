# MeltR

Automated fitting of RNA/DNA absorbance melting curves and fluorescence binding isotherms in R

Jacob P. Sieg\textsuperscript{[1],[2]}, Philip C. Bevilacqua\textsuperscript{[1],[2],[3]}

\textsuperscript{[1]}Department of Chemistry, Pennsylvania State University, University Park, PA 16802.

\textsuperscript{[2]}Center for RNA Molecular Biology, Pennsylvania State University, University Park, PA 16802.

\textsuperscript{[3]}Department of Biochemistry and Molecular Biology, Pennsylvania State University, University Park, PA 16802.

Last update, 12/02/2022

# 1 Overview

MeltR is a R package, written by Jacob Sieg, that fits nucleic acid folding data to molecular models to obtain thermodynamic parameters. MeltR automates the trivial but time-consuming tasks associated with fitting nucleic acids thermodenaturation data, leading to facile conversion of raw data into useful thermodynamic parameters.

MeltR was inspired by Meltwin.\textsuperscript{[1]} MeltR and Meltwin have the same utility: easy and consistent fitting to obtain thermodynamic parameters. The main drawback of MeltR is that it is ran from your R console, whereas Meltwin has a graphical-user-interface. However, the MeltR syntax is not complicated, and MeltR has other advantages: (1) A current versions of MeltR can be downloaded from GitHub by entering two lines of code in your R console, whereas Meltwin has been out of support for years. (2) MeltR supports fitting fluorescence binding isotherms to obtain thermodynamic parameters. (3) MeltR can be ran in bulk. (4) MeltR can be run on any operating system that is compatible with R. (5) Anecdotally, MeltR is more robust than Meltwin and requires less input from the user.

The core of MeltR is the “meltR.A” function for fitting absorbance melting curves and the “meltR.F” function for fitting fluorescence binding isotherms.

The creators of MeltR are interested in improving MeltR by adding new molecular models. Areas of interest include multistate models and models to describe G-quadruplex RNA. Please email Jacob Sieg at jus841@psu.edu for suggestions.  

# 2 MeltR installation

MeltR is written in R. You can find instructions for installing R here: "https://www.r-project.org/".

MeltR can be installed from "https://github.com/JPSieg/MeltR" by typing two lines into your R console.

```{r}
install.packages("devtools")
devtools::install_github("JPSieg/MeltR")
```

MeltR is written entirely in base R and requires no other dependencies. However, Rstudio, the tidyverse, and ggrepel, are recommended software and packages for running MeltR. You can find instructions for installing Rstudio here: "https://www.rstudio.com/products/rstudio/download/". The free version is excellent. Tidyverse and ggrepel can be installed from your R console.

```{r}
install.packages("tidyverse")
install.packages("ggrepel")
```

# 3 Theory

## 3.1 Fitting fluorescence isotherms

MeltR can obtains thermodynamic parameters from fluorescence binding isotherms for heteroduplex DNA and RNA using the non-linear regression function in base R, "nls". In this strategy, a quencher labeled strand is titrated into different wells in a qPCR plate containing a constant concentration of fluorophore labeled strand. The fluorophore labeled strand binds to the quencher labeled strand, reducing the fluorescence emission (E) and resulting in an apparent fluorescence binding isotherm, where the shape of the curve is determined by Equation 1.

\begin{equation} \label{eqn}
E = F_{max} + (F_{min} - F_{max})*F(K_{D}, [A]_{T}, [B]_{T})
	\end{equation}

Where F\textsubscript{max} is the fluorescence emission of unbound fluorophore labeled strand (A), F\textsubscript{min} is the fluorescence emission of the fluorophore labeled strand completely bound to a quencher labeled strand (B), and f(K\textsubscript{D}, [A]\textsubscript{T}, [B]\textsubscript{T}) is the mole fraction of the bound fluorophore labeled strand (AB) to the total fluorophore labeled strand.

\begin{equation} \label{eqn}
f(K_{D}, [A]_{T}, [B]_{T}) = \frac{[AB]}{[A]_{T}}
	\end{equation}

f(K\textsubscript{D}, [A]\textsubscript{T}, [B]\textsubscript{T})  is a function of the total fluorophore labeled strand ([A]\textsubscript{T}) and the total quenceher labeled strand ([B]\textsubscript{T}), and the K\textsubscript{D} is given by the the expression:

\begin{equation} \label{eqn}
K_{D} = \frac{[A][B]}{[AB]}
	\end{equation}

F(K\textsubscript{D}, [A]\textsubscript{T}, [B]\textsubscript{T}) can be determined by solving the K\textsubscript{D} expression to obtain:


\begin{equation} \label{eqn}
F(K_{D}, [A]_{T}, [B]_{T})  = \frac{(K_{D}+[A]_{T}+[B]_{T}) - \sqrt{{(K_{D}+[A]_{T}+B]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}
	\end{equation}

Thus, the K\textsubscript{D} at each temperature was determined by fitting isotherms at each temperature to Equation 5, obtained by plugging Equation 4 into Equation 1.

\begin{equation} \label{eqn}
E = F_{max} + (F_{min} - F_{max})\frac{(K_{D}+[A]_{T}+[B]_{T}) - \sqrt{{(K_{D}+[A]_{T}+B]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}
	\end{equation}

Thermodynamic parameters for helix formation were extracted by the Van't Hoff relationship (Equation 6) between the K\textsubscript{D} and the temperature using the two methods described below.

\begin{equation} \label{eqn}
ln(K_{D}) = \frac{dS}{R} - \frac{dH}{RT} 
	\end{equation}

The dS is the entropy of helix formation, the dH is the enthalpy of helix formation, R is the gas constant, and T is the temperature in Kelvin. For a single experiment, all samples should be diluted from the same fluorophore and quencher stocks, so that MeltR can deal with systematic uncertainty in strand concentrations. 

### 3.1.1 Fluorescence data preprocessing

MeltR performs the following data preprocessing steps before fitting:

1.) The fluorophore labeled strand concentration is optimized using the concentration optimization algorithm, described in more detail in section 3.1.3.

2.) The T\textsubscript{m} of each sample is estimated using the first derivative of fluorescence emission as a function of temperature. First derivatives are calculated using polynomial regression, where the data are approximated with a 20th order polynomial. Then, the analytical first derivative of the polynomial was determined using calculus. The approximate T\textsubscript{m} (where 50% of the nucleic acid is single stranded) is calculated by finding the maximum of the first derivative curve to a precision of less than 0.1 <span>&#176;</span>C. This T\textsubscript{m} is unreliable because fluorescence baselines vary unpredictable with temperature. However, the T\textsubscript{m} is useful for qualitative comparison of stability between conditions.  

2.) Isotherms are fit to Equation 5 to determine K\textsubscript{D} and error in the K\textsubscript{D} at each temperature using "nls" in base R. Initial values for F\textsubscript{max} and F\textsubscript{min} are estimated by taking the mean of the 20% highest readings in each isotherm and the 20% lowest readings respectively. Initial values for the K\textsubscript{D} are provided by the user, by default 0.1 nM. 

3.) K\textsubscript{D}s are filtered by magnitude and error, according to user specifications, to determine which isotherms are most reliable. First, K\textsubscript{D}s outside of a user specified range (10 to 1000 nM by default) are thrown out. Second, K\textsubscript{D}s are ranked by the error in the K\textsubscript{D}. K\textsubscript{D}s that are below a user specified error quantile are thrown out. The default K\textsubscript{D} error quantile file is 0.25, meaning the algorithm will keep the to 25% most accurate K\textsubscript{D}s, after filtering by magnitude.  

4.) The most reliable K\textsubscript{D}s and temperatures are passed to Method 1 to make Van't Hoff plots and to calculate helix formation energies.

5.) Fluorescence data from the most reliable K\textsubscript{D}s and temperatures are passed to Method 2 for global fitting and to calculate helix formation energies. 

### 3.1.2 Concentration optimization algorithm.

We have determined that the accuracy of fit results is highly dependent on errors in the predicted concentration of fluorophore and quencher RNA strands in stock solutions. Fit accuracy is dependent on the mole ratio of fluorophore and quencher labeled RNA in stock solutions, but not dependent on the total magnitude of both fluorophore and quencher labeled RNA concentrations in the stocks. Thus, the concentration optimization algorithm in MeltR does not need to find exact concentrations to generate accurate helix formation energies, just the mole ratio (X) of fluorophore and quencher labeled RNA in samples that are predicted to have equal concentrations of fluorophore and quencher, for example predicted 200 nM FAM-RNA and 200 nM RNA-BHQ1.

MeltR uses a fluorescence binding isotherm from a temperature where the K\textsubscript{D} is more than 10 times less than the fluorophore labeled strand concentration. Under these conditions, the shape of the curve will be independent of the K\textsubscript{D} , and the isotherm can be used as a Job plot. For example, at a 200 nM predicted fluorophore labeled strand concentration, the shape of the binding curve will be independent of K\textsubscript{D}  if the K\textsubscript{D}  is less than 10 nM. The curve will resemble a hockey stick (Figure 1A) composed of two straight lines. The first line, where [A]\textsubscript{T} > [B]\textsubscript{T}, will decrease as the [B]\textsubscript{T} is increased and the fluorophore labeled strand is saturated. The second line, where [A]\textsubscript{T} < [B]T\textsubscript{T}, will be horizontal as the [B]\textsubscript{T} is increased because the fluorophore labeled strand is saturated. The intersection of the first and second line will occur at:

\begin{equation} \label{eqn}
[A]_{T} = [B]_{T}
	\end{equation}

![Job plot used by the concentration algorithm. Data are modeled with a KD of 0.1 nM and a FAM-RNA strand concentration of 250 nM. (A) Isotherms resemble a job plot, where the intersection of Line 1 and Line 2 can be used to determine X. (B) The same modeled data with the concentration algorithm starting at X = 1 and ending at X = 0.8.](https://user-images.githubusercontent.com/63312483/166508109-8337d818-6eed-4b08-a33c-92794a14a95a.svg)

The absolute concentration of [A]\textsubscript{T} and [B]\textsubscript{T} cannot be known with precicion. However, the mole ratio of the error (X) in the fluorophore and quencher stocks can be estimated from this data point at Equation 8.

\begin{equation} \label{eqn}
\frac{[A]_{T-estimated}}{X} = [B]_{T}
	\end{equation}

Where [A]\textsubscript{T-estimated} is the estimated total A concentration, or the predicted concentration based on of the stock, and X is given by Equation 9 based on the actual concentration determined from the experiment.

\begin{equation} \label{eqn}
X = \frac{[A]_{T}}{[B]_{T}}
	\end{equation}

MeltR fits an overdetermined isotherm to a binding curve (selected by the user but by default the isotherm collected at the lowest temperature) to a modified version of Equation 5 to determine X (Figure 1 B).

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = F_{max} + (F_{min} - F_{max})*\frac{(K_{D}+[A]_{T-estimated}+[B]_{T}X) - \sqrt{{(K_{D}+[A]_{T-estimated}+B]_{T}X)}^2 - 4[A]_{T-estimated}[B]_{T}X}}{2[A]_{T-estimated}}$}
	\end{equation}

By default, the fit to determine X allows X and K\textsubscript{D} to float. The user should also use an argument called "low_K", to set the Kd in the optimization fit to several K\textsubscript{D} that are more than 10 times less than the fluorophore labeled strand concentrations, as described in Section 4.1.3. The user should then inspect X from several iterations of the optimization algorithm set to different "low_K" values to make sure it is similar to the iteration that allows the K\textsubscript{D} to float.

[A]\textsubscript{T} is then corrected with Equation 11.

\begin{equation} \label{eqn}
[A]_{T}  = \frac{[A]_{T-estimated}}{X}
	\end{equation}

### 3.1.3 Method 1 Van't Hoff plot

Method 1 fits K\textsubscript{D}s that were passed from the preprocessing steps to Equation 6. Initial estimates for the dH and dS of helix formation are provided by the user in kcal as a list using the "vh_start" argument. The default is -70 kcal/mol and -0.180 kcal/mol/K for the enthalpy and entropy respectively, and will work for most helices. The free energy of helix formation at 37 <span>&#176;</span>C (dG) was calculated from the dH and dS values provided by the fit (Equation 12).

\begin{equation} \label{eqn}
dG = dH - 310.15*dS \qquad (12)
	\end{equation}

The standard error of the dG was propagated from the standard error of dH and the dS from the fit, and the covariation between the two values.

\begin{equation} \label{eqn}
SE_{dG} = \sqrt{ {SE_{dH}}^2 + {(310.15* SE_{dS})}^2 + 2* 310.15* Covar_{dH,dS}}
	\end{equation}

Helix formation energies are traditionally reported in terms of the association constant.

\begin{equation} \label{eqn}
K = \frac{[AB]}{[A][B]} = \frac{1}{K_{D}}
	\end{equation}

The dHs and dSs reported by MeltR are also reported in terms of the association constant, which is obtained by multiplying helix formation energies from fitting K\textsubscript{D} by negative one. 

### 3.1.4 Method 2 Global fit

Method 2 globally fits fluorescence data that were passed from the preprocessing to Equation 15, obtained by plugging Equation 6 into Equation 5. F\textsubscript{max} and F\textsubscript{min} is allowed to float between temperatures and dH and dS are fixed between temperatures. Starting values for F\textsubscript{max}, F\textsubscript{min}, dH, and dS are obtained from the results for individual fits. 

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = F_{max} + (F_{min} - F_{max})*\frac{(\exp{(\frac{dS}{R} - \frac{dH}{RT})}+[A]_{T}+[B]_{T}) - \sqrt{{(\exp{(\frac{dS}{R} - \frac{dH}{RT})}+[A]_{T}+B]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}$}
	\end{equation}

The dG and the error in the dG were calculated as in method 1.

## 3.2 Fitting abosorbance melting curves

### 3.2.1 Absorbance data preprocessing

MeltR can obtain thermodynamic parameters from absorbance melting curves for heteroduplex, homoduplex (selfcomplementary), and monomolecular self-structured DNA and RNA using the method pioneered by Tinoco and colleagues and refined by Turner and colleagues.\textsuperscript{[1-3]}

MeltR performs the following data preprocessing steps before fitting:

1.) RNA and DNA extinction coefficients are calculated from the sequence using the extinction coefficient data from Tinoco and colleagues.\textsuperscript{[2]}

2.) The background absorbance of a user specified blank is scaled to the pathlength of the cuvette and subtracted from each curve.

3.) The total strand concentration (C\textsubscript{t}) is calculated at a user specified temperature (Default = 90 <span>&#176;</span>C). For homoduplexes and monomolecular-self structured RNA, where there is only one nucleic acid in the sample, the C\textsubscript{t} is the total concentration of the single strand. For heteroduplexes, the C\textsubscript{t} is the total concentration of strand 1 plus the total concentration of strand 2.

4.) High and low temperature data are trimmed (to the users specification) to ensure linear baselines.

5.) First and second derivatives are taken using polynomial regression. First, the data are approximated using a 20th order polynomial. Then, the first and second derivatives of the polynomial are determined analytically using calculus. The approximate T\textsubscript{0.5} (the approximate melting temperature, T\textsubscript{m}, where 50% of the nucleic acid is single stranded) is calculated by finding the maximum of the first derivative curve and the T\textsubscript{0.75} (the approximate temperature where 75% of the nucleic acid is single stranded) is calculated by finding the minimum of the second derivative (Figure 2), to a precision of less than 0.1 <span>&#176;</span>C.

![Absorbance melting curve derivative determined with polynomial regression.](https://user-images.githubusercontent.com/63312483/166296578-341d3aac-05a9-4176-9a72-0831dc0487cf.svg)

6.) Initial parameter estimates are calculated for each curve. The initial values for slopes and intercepts of the baselines are estimated by fitting absorbance values that are greater than the 75th quantile for the uper baselines, and fitting aborbance values that are lower than the 25th quantile for the lower baseline, to y = mx + b. Initial values for the enthalpy are determined using the T\textsubscript{0.5} and T\textsubscript{0.75} (in Kelvin) from first and second derivative curves. MeltR uses equation 16, equation 17, and equation 18, for heteroduplex, homoduplex, and monomolecular self-structured DNA and RNA melting curves.   

\begin{equation} \label{eqn}
dH = -0.007*(\frac{1}{T_{0.5}} - \frac{1}{T_{0.75}})
	\end{equation}

\begin{equation} \label{eqn}
dH = -0.0044*(\frac{1}{T_{0.5}} - \frac{1}{T_{0.75}})
	\end{equation}

\begin{equation} \label{eqn}
dH = -0.0032*(\frac{1}{T_{0.5}} - \frac{1}{T_{0.75}})
	\end{equation}

The initial T\textsubscript{m} is estimated from the T\textsubscript{0.5}, determined in step 5.

### 3.2.2 Thermodynamic models for duplex formation

Thermo-dynamic parameters for helix formation are obtained using a Van't Hoff model:

\begin{equation} \label{eqn}
lnK = \frac{dS}{R} - \frac{dH}{RT}
	\end{equation}

where dS is the entropy change, dH is the enthalpy change, R is the gas constant in kcal/mol, T is the temperature in Kelvin, and K is the equillibrium constant given by Equation 20, Equation 21, and Equation 22 for heteroduplexes, homoduplexes, and monomolecular self-structured RNA respectively.

\begin{equation} \label{eqn}
K = \frac{[AB]}{[A][B]}
	\end{equation}

\begin{equation} \label{eqn}
K = \frac{[AA]}{[A]^2}
	\end{equation}

\begin{equation} \label{eqn}
K = \frac{[F]}{[U]}
	\end{equation}

For Equation 20 and Equation 21, [A] and [B] are the concentration of different strands, [AB] is the concentration of strand A in a duplex with strand B, and [AA] is the concentration of a self complementary strand A in a duplex with another self complementary strand A. For Equation 22, [F] is the concentration of a monomolecular self-structured RNA in the folded state and [U] is the concentration of a monomolecular self-structured RNA in the unfolded state.

MeltR uses three methods based on the Van't Hoff equation to calculate thermodynamic parameters: (1) fitting melting curves individually, (2) fitting the thermodenaturation point as a function of temperature, and (3) Global fitting melting curves.

### 3.2.3 Method 1 fitting melting curves individually

Method 1 fits the absorption as a function of temperature for each sample individually. Base lines are modeled as y= mx + b for the absorbance of the unfolded-single stranded and folded-duplex states (Equation 23). 

\begin{equation} \label{eqn}
A = mT + b
	\end{equation}

The absorption of each sample as a function of temperature is also a function of the fraction of RNA in the folded-duplex state (DS), as a function of temperature f(T), given by Equation 24. SS represents the RNA in the single-stranded state.

\begin{equation} \label{eqn}
A = (m_{DS}T + b_{DS})f(T) + (m_{SS}T + b_{SS})(1-f(T))
	\end{equation}

f(T) is variable, calculated by the analytic solution of the binding constant. MeltR uses Equation 25 for heteroduplexes, Equation 26 for homoduplexes, and Equation 27 monomolecular self-structured RNA.

\begin{equation} \label{eqn}
f(T) = \frac{\frac{2}{K(T)*Ct} + 2 - \sqrt{(\frac{2}{K(T)*Ct} + 2)^2 - 4}}{2}
	\end{equation}

\begin{equation} \label{eqn}
f(T) = \frac{\frac{1}{2*K(T)*Ct} + 2 - \sqrt{(\frac{1}{2*K(T)*Ct} + 2)^2 - 4}}{2}
	\end{equation}

\begin{equation} \label{eqn}
f(T) = \frac{K(T)}{1 + K(T)}
	\end{equation}

Where C\textsubscript{t} is the total strand concentration. K(T) is the equilibrium constant as a function of temperature, given by Equation 28 for heteroduplexes, Equation 29 for homoduplexes, and Equation 30 for monomolecular self-structured RNA.

\begin{equation} \label{eqn}
K(T) = \exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{4}{Ct}))}
	\end{equation}

\begin{equation} \label{eqn}
K(T) = \exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{1}{Ct}))}
	\end{equation}

\begin{equation} \label{eqn}
K(T) = \exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}))}
	\end{equation}

Note, K(T) is in terms of T\textsubscript{m} and C\textsubscript{t}, instead of the dS, to increase the ease of estimating initial parameters for non-linear regression and to increase the robustness of the nls algorithm. Briefly, the dS is replaced with T\textsubscript{m} and C\textsubscript{t} by solving Equation 19 for the dS at the T\textsubscript{m}, where f(T) = 0.5.

Thus, method 1 fits absorbtion versus temperature for each sample to equations 31, 32, and 33 to determine thermodynamic prameters for heteroduplexes, homoduplexes, and monomolecular self-structured RNA respectively.

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = (m_{DS}T + b_{DS})\frac{\frac{2}{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{4}{Ct}))}*Ct} + 2 - \sqrt{(\frac{2}{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{4}{Ct}))}*Ct} + 2)^2 - 4}}{2} + (m_{SS}T + b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{4}{Ct}))}*Ct} + 2 - \sqrt{(\frac{2}{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{4}{Ct}))}*Ct} + 2)^2 - 4}}{2})$}
	\end{equation}

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = (m_{DS}T + b_{DS})\frac{\frac{1}{2*\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{1}{Ct}))}*Ct} + 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{1}{Ct}))}*Ct} + 2)^2 - 4}}{2} + (m_{SS}T + b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{1}{Ct}))}*Ct} + 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}) + ln(\frac{1}{Ct}))}*Ct} + 2)^2 - 4}}{2})$}
	\end{equation}

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = (m_{DS}T + b_{DS})\frac{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}))}}{1 + \exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}))}} + (m_{SS}T + b_{SS})(1-\frac{\exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}))}}{1 + \exp{(\frac{dH}{R}(\frac{1}{T_{m}} - \frac{1}{T}))}})$}
	\end{equation}

Free energy at 37 <span>&#176;</span>C (dG) is calculated from the dH and entropy (dS) of helix formation 

\begin{equation} \label{eqn}
dG = dH - 310.15*dS
	\end{equation}

The dS of helix formation is calculated from the dH and the Tm.

For heteroduplexes:

\begin{equation} \label{eqn}
dS = \frac{dH}{Tm} + Rln(\frac{4}{Ct})
	\end{equation}

For homoduplexes:

\begin{equation} \label{eqn}
dS = \frac{dH}{Tm} + Rln(\frac{1}{Ct})
	\end{equation}

For monomolecular self-structured RNA/DNA:

\begin{equation} \label{eqn}
dS = \frac{dH}{Tm}
	\end{equation}

The mean and the standard deviations are reported for thermodynamic parameters.

### 3.2.4 Method 2 fitting the T\textsubscript{m} as a function of C\textsubscript{t}

Method 2 uses nls to fit the relationship between 1/T\textsubscript{m} and C\textsubscript{t} to Equation 38 and 39 for heteroduplexes and homoduplexes, respectively. The T\textsubscript{m} of monomolecular, self-structured RNA is independent of C\textsubscript{t}, so Method 2 cannot be used, and MeltR automatically turns Method 2 off.

\begin{equation} \label{eqn}
\frac{1}{T_{m}} =  \frac{R}{dH}*lnC_{t} + \frac{dS}{dH} - R*ln(4)
	\end{equation}

\begin{equation} \label{eqn}
\frac{1}{T_{m}} =  \frac{R}{dH}*lnC_{t} + \frac{dS}{dH}
	\end{equation}
	
The T\textsubscript{m} for each sample calculated from method 1 is used for method 2 by default. The propagated error in the T\textsubscript{m} from the method 1 fits can be used as weights, propagated from the method 1 fits using Equation 40. The weighted
regression can be turned turned of by the user but the default setting is to perform a non-weighted regression for method 2.

\begin{equation} \label{eqn}
Weight = \frac{{T_{m}}^2}{SE_{T_{m}}}
	\end{equation}
	
MeltR can use two other T\textsubscript{m} determination methods and neither can use a weighted regression. The first uses the linear baseline parameters from Method 1 to calculate the f(T) using Equation 41.

\begin{equation} \label{eqn}
 f(T) = \frac{A -(m_{SS}T + b_{SS})}{(m_{DS}T + b_{DS}) + (m_{SS}T + b_{SS})}
	\end{equation}

F(T) is approimatelty linear in the range of 0.4 to 0.6. Thus, F(T in {0.4 to 0.6}) was fit with y = mT + b, and solved using y = 0.5 to accurately determine the melting temperature for each C\textsubscript{t}. Alternatively, Method 2 can use the T\textsubscript{m} estimated using polynomial regression followed by the first derivative analysis (step 5 of data processing)


Free energy at 37 <span>&#176;</span>C (dG) is calculated from the dH and entropy (dS) of helix formation directly from the fit.

\begin{equation} \label{eqn}
dG = dH - 310.15*dS
	\end{equation}

Error in the dG is calculated by propagating error in the fit terms dH and dS.

\begin{equation} \label{eqn}
SE_{dG} = \sqrt{ {SE_{dH}}^2 + {(310.15* SE_{dS})}^2 - 2*310.15*Covar_{dH, dS}}
	\end{equation}

### 3.2.5 Method 3 Global fitting

Method 3 fits all curves to Equations 44, 45, and 46, simultaneously in a global fit, for heteroduplexes, homoduplexes, and mono-molecular self-structured RNA/DNA respectively. In this global fit, Equations 31, 32, and 33, are rearranged to be in terms of the dS instead of the Tm. 

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = (m_{DS}T + b_{DS})\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2)^2 - 4}}{2} + (m_{SS}T + b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2)^2 - 4}}{2})$}
	\end{equation}

\begin{equation} \label{eqn}
\resizebox{\hsize}{!}{$E = (m_{DS}T + b_{DS})\frac{\frac{1}{2\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2)^2 - 4}}{2} + (m_{SS}T + b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} + 2)^2 - 4}}{2})$}
	\end{equation}

\begin{equation} \label{eqn}
E = (m_{DS}T + b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 +\exp{(\frac{dS}{R} - \frac{dH}{RT})}} + (m_{SS}T + b_{SS})(1-\frac{\exp{(\frac{dS}{R} - \frac{dH}{RT})}}{1 + \exp{(\frac{dS}{R} - \frac{dH}{RT})}})
	\end{equation}

The baselines are allowed to vary but dHs and dSs are constrained to a single value for all curves. For global fitting, the slopes and intercepts of the fits from Method 1 are used as initial parameter estimates for the slopes and intercepts of the global fit, and the average of the dHs and dSs from Method 1 are used as initial parameter estimates for the dH and dS.

The dG and error in the dG is calculated using the same equations as Method 2.

### 3.2.6 T\textsubscript{m} at 0.1 mM

The T\textsubscript{m} at 0.1 mM is not a true thermodynamic parameter and included in the output for historical reasons. It is an expectation maximized value, meaning it is estimated from a fit. The T\textsubscript{m}at 0.1 mM is calculated using Equation 47, obtained by rearranging equation 35, 36, and 37.

\begin{equation} \label{eqn}
{T_{m}}^{0.1 mM} = \frac{dH}{dS - R*lnB}
\end{equation}

B is 4/C\textsubscript{t}, 1/C\textsubscript{t}, or 1 for heteroduplexes, homoduplexes, and monomolecular self structured RNA, respectively.

The propagation of error to the T\textsubscript{m} at 0.1 mM was calculated with Equation 48.

\begin{equation} \label{eqn}
SE_{{T_{m}}^{0.1 mM}} = \sqrt{ {(\frac{SE_{dH}}{dS - RlnB})}^2 + {(\frac{dH* SE_{dS}}{{(dS - RlnB)}^2})}^2 - 2\frac{dH*Covar_{dH,dS}}{{(dS - RlnB)}^3}}
\end{equation}

Covariance is estimated using Equation 49 for method 1, where the u\textsubscript{dH} and u\textsubscript{dS} represent the mean of the individual fits. Covariance is extracted from the fit for method 2 and 3.

\begin{equation} \label{eqn}
Covar_{dH,dS} = \frac{1}{N}\sum{(dH_{i} - u_{dH})(dH_{i} - u_{dH})}
\end{equation}

### 3.2.7 Autobaseline trimming using the BLTrimmer

MeltR provides an automated baseline trimming function, “BLtrimmer”, which works on three principles:

1. Many combinations of trimmed baselines are generated and fit.
2. The best baseline combinations produce the most internally consistent folding thermodynamic parameter estimates.
3. Thermodynamic parameter from the best combinations should be treated as an ensemble of equally feasible estimates.

Thus, the BLTrimmer function explores the dependence of fit parameters on baseline trimming, starting from an existing meltR.A fit object. The BLTrimmer uses three steps: Rapid testing of baseline combinations, assessment, and a slower treatment of an ensemble of baselines that produce internally consistent thermodynamic energies.

Baseline testing: For the first step, the BLtrimmer generates and then fits a large number of baseline combinations. The user can generate baseline combinations using a fixed or a floating baseline trim. The fixed method uses the same baseline length for each sample to test how baseline length effects the subsequent thermodynamic parameters. The floating method allows different samples to have different baseline lengths. The floating method generates a large number of possible baseline combinations. Thus, the user can only test a small subset of these baseline combinations using a reasonable amount of computational time. Thus, a randomly selected subset of baseline combinations are tested using a fast, lightweight fit. The BLTrimmer fits each sample using Method 1, using the global fit from the meltR.A object as initial guesses for the non-linear regression step. Using these good initial guesses reduces the time the nls algorithm needs to find a solution. Next, the BLTrimmer fits each sample using Method 2. dH values from each method and baseline combination are then passed to the assessment step.

Baseline assessment: The user can choose 1 of 3 assessment methods to identify the baseline combinations that produce the most internally consistent thermodynamic parameters.

Assessment Method 1 tests how well the dH values from each sample agree with each other. The normalized standard deviation SE\textsubscript{dH1}' is calculated for each baseline combination, where SE\textsubscript{dH1} and u\textsubscript{dH1} are the standard deviation and mean dH, respectively from the samples in Method 1.

\begin{equation} \label{eqn}
SE_{dH1}' =  \frac{|SE_{dH1}|}{|u_{dH1}|}
	\end{equation}

SE\textsubscript{dH1}' is then ranked into quantiles on a scale from 0 to 1. For example, if a baseline combination has the 100th smallest SE\textsubscript{dH1}' in a set of 1000 baseline combinations, its quantile (Q\textsubscript{dH1}) is 0.1. Then, baseline combinations with Q\textsubscript{dH1} values smaller than a user specified value are passed to the ensemble analysis.

Assessment Method 2 tests how well the dH values from Method 1 and 2 to agree using equation 49, where dH2 is produced by Method 2.

\begin{equation} \label{eqn}
SE_{dH1,2}' =  2\frac{|SE_{dH1} - dH2|}{|u_{dH1} + dH2|}
	\end{equation}

SE\textsubscript{dH1,2}' is then ranked into quantiles on a scale from 0 to 1 to produce Q\textsubscript{dH1,2} values and baseline combinations. Q\textsubscript{dH1,2} values smaller than a user specified value are passed to the ensemble analysis.

Assessment Method 3 combines assessment methods 1 and 2 using a statistic called the error distance (Q\textsubscript{D}).

\begin{equation} \label{eqn}
SE_{D} = \sqrt{{Q_{dH1}}^2 - {Q_{dH1,2}}^2}
	\end{equation}

SE\textsubscript{D} (equation 50) is then ranked into quantiles on a scale from 0 to 1 to produce  Q\textsubscript{D} values and baseline combinations. Q\textsubscript{D} values smaller than a user specified value are passed to the ensemble analysis.

Ensemble analysis: Baselines that pass the assessment criterion are treated as an ensemble of equally feasible baseline combinations. Each baseline combination is passed to meltR.A and fit, which takes considerably longer than the lightweight analysis in the testing step because initial guesses are predicted for each trim using the protocol described for meltR.A. The average of resulting thermodynamic parameters are reported with 95% confidence intervals.

# 4 Running MeltR

In this section, we cover how to use MeltR in your R console. If you have not already, read section 3 to understand the theory underlying the results of MeltR. Section 4 will cover MeltR usage and how to avoid pitfalls. The most common error with MeltR is a user attempting to fit data that is inconsistent with the underlying model, either a fluorescence isotherm or a absorbance melting curve with a non-standard shape or specifying an incorrect molecular model. In the case of data with a non-standard shape, the aberrant data will need to be filtered out data set prior to fitting the set with MeltR. While MeltR has no dependencies other than base R 4.1.3, data wrangling and plotting functions in the "tidyverse" package\textsuperscript{[4]} are highly recommended, along with the "ggrepel" package\textsuperscript{[5]}, for data quality checks and filtering. To begin, open a new R session in the proper directory. Load relevant packages into your memory.

```{r}
library(MeltR)
library(tidyverse)
library(ggrepel)
```

Help documentation for MeltR functions can be pulled up using standard R commands.

```{r}
?meltR.F
?meltR.A
```

## 4.1 Fitting fluorescence binding isotherms

### 4.1.1 Formatting fluorescence data for MeltR

Data should be formatted into a comma separated value (“.csv”) text file with carefully labeled columns (Figure 3). There should be a “Well” column where numbers or character strings specify what well in a microplate the data came from, a “Reading” column that specifies the reading a data point comes from, a “Temperature” column that specifies the temperature where the data was recorded, a “B” column which specifies an approximate quencher labeled strand concentration in nM, an “A” column which specifies an approximate fluorophore labeled strand concentration in nM, and an “Emission” column containing the fluorescence emission intensity. 

![Formatting fluorescence binding isotherm datat in a csv file for MeltR with Excel. B and A are the concentration of quencher and fluorophore labeled RNA strands, respectively, in nM. Temperature is in <span>&#176;</span>C.](https://user-images.githubusercontent.com/63312483/165355843-40b03fbc-7d89-429e-8252-1b65ea6fced1.png)

Two common pitfalls occur when formatting data for MeltR, usually in Excel. The first is incorrect column names, even incorporation of an extra space, so that MeltR cannot recognize relevant data when it is read into R. The second is incorporation of characters into data columns when values are missing. If a data point is missing, leave the cell blank. Do not write something like "NA".

### 4.1.2 Reading data into a R data frame

Comma separated value (“.csv”) data can be read into a R data frame for MeltR using the “read.csv” function that is included in base R.

```{r}
df = read.csv("path/file_name.csv")
```

Data can be checked with the "View" function.

```{r}
View(df)
```

Note, the columns in “df” must be named correctly. If the columns are not named correctly, MeltR cannot recognize them. You can rename the columns of a data frame using:

```{r}
colnames(df) = c("Helix", "Well", "Reading", "Temperature", "B", "A", "Emission")
```

### 4.1.3 Applying “meltR.F” to obtain thermodynamic parameters

For this tutorial, we will use sample data included in MeltR. First check that MeltR is in your memory and find out what is in the sample data.

```{r}
?df.fluor.data
df = df.fluor.data
```

In general, it is a good idea to plot your data before you start fitting. First, I want to check the quality of low temperature isotherms. Using the tidyverse and ggrepel.

```{r}
head(df)
ggplot(df %>% filter(Reading == 1), aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()
```

The result is a really nice binding isotherm (Figure 4). One should check a few more readings by changing the value in the filter command, for example 20, 60, 80, etc... 

![Fluorescence isotherm example from data included in MeltR. Labels represent the well in a 96 well plate.](https://user-images.githubusercontent.com/63312483/165367119-ed4e15df-4c08-46d7-82c6-03b386f61395.svg)

If you observe one or two obvious outliers in the data set, it is reasonable to remove them using the filter command. However, this data is excellent and needs no filtering.

```{r}
df %>% filter(!Well %in% c("A1", "C3"))
```

If you are satisfied with the the data, we can move on to fitting. I will first inspect the help file for meltR.F.

```{r}
?meltR.F
```

The only argument that needs set is the data frame.

```{r}
meltR.F(df)
```

You will see the following output in your console.

```{r}
[1] "Van't Hoff"
[1] "accurate Ks = 16"
        Method         H     SE.H         S      SE.S         G       SE.G
1    1 VH plot -48.37473 1.048500 -121.6157  3.357021 -10.65562 0.01033199
2 2 Global fit -48.88546 6.302188 -123.2438 20.157522 -10.66140 0.06564137
    K_error         R   Kd.opt
1 0.2916302 0.7797618 6.165553
2 0.2916302 0.7797618 6.165553
[1] "Fractional error between methods"
           H        S           G
1 0.01050236 0.013298 0.000542719
```

Note, the concentration optimization algorithm, by default allowing the Kd.opt to float, identified an optimal R of ~0.78. Since the FAM-RNA strand in each well is estimated at 200 nM, the concentration optimization algorithm will adjust the concentration to 200/0.78 = 256 nM. We will test the robustness of this estimate by constraining the Kd.opt to a range of K\textsubscript{D}s that are more than 10 times less than the FAM-RNA concentration, using the "low_K" argument. For example, 1, 0.5, 0.1, and 0.05 nM. Note, the mole ratio "R" was labeled "X" in the theory section to avoid confusion with the gas constant.

```{r}
> meltR.F(df, low_K = 1)
[1] "Van't Hoff"
[1] "accurate Ks = 11"
        Method         H       SE.H         S      SE.S         G        SE.G
1    1 VH plot -60.55018  0.7291223 -160.1005  2.322061 -10.89501 0.009578498
2 2 Global fit -60.70845 12.4962966 -160.6024 39.782555 -10.89761 0.167951340
    K_error         R Kd.opt
1 0.3333168 0.7231654      1
2 0.3333168 0.7231654      1
[1] "Fractional error between methods"
            H           S            G
1 0.002610529 0.003130137 0.0002389311
> meltR.F(df, low_K = 0.5)
[1] "Van't Hoff"
[1] "Van't Hoff"
[1] "accurate Ks = 11"
        Method         H       SE.H         S      SE.S         G        SE.G
1    1 VH plot -61.48056  0.7017499 -163.0171  2.234887 -10.92082 0.009218896
2 2 Global fit -61.64895 12.7501200 -163.5512 40.589141 -10.92354 0.171775450
    K_error         R Kd.opt
1 0.3374478 0.7153701    0.5
2 0.3374478 0.7153701    0.5
[1] "Fractional error between methods"
           H           S            G
1 0.00273514 0.003271292 0.0002491787
> meltR.F(df, low_K = 0.1)
[1] "Van't Hoff"
[1] "accurate Ks = 14"
[1] "Van't Hoff"
[1] "accurate Ks = 10"
        Method         H       SE.H         S      SE.S         G        SE.G
1    1 VH plot -63.04998  0.6809662 -167.9641  2.167052 -10.95592 0.009325515
2 2 Global fit -63.18785 14.8873698 -168.4012 47.361747 -10.95822 0.207758608
    K_error         R Kd.opt
1 0.3446286 0.7080886    0.1
2 0.3446286 0.7080886    0.1
[1] "Fractional error between methods"
            H           S            G
1 0.002184357 0.002599086 0.0002100154
> meltR.F(df, low_K = 0.05)
        Method         H      SE.H         S      SE.S         G        SE.G
1    1 VH plot -63.17200  0.677299 -168.3465  2.155382 -10.95934 0.009275292
2 2 Global fit -63.31084 14.925896 -168.7867 47.484124 -10.96165 0.208349999
    K_error         R Kd.opt
1 0.3451437 0.7070892   0.05
2 0.3451437 0.7070892   0.05
[1] "Fractional error between methods"
            H           S            G
1 0.002195498 0.002611534 0.0002110328
```

The optimization algorithm, using a constrained K\textsubscript{D}, adjusts the FAM-RNA concentration to at most 281, about 10% difference from the optimization algorithm using a floating K\textsubscript{D}. Thus, the default concentration optimization is robust. Note that there is considerable variance in the enthalpies, entropies, and free energies from the various constrained K\textsubscript{D}s and the unconstrained fit (about 20 kcal/mole in terms of the enthalpy). We will demonstrate how to refine these thermodynamic parameters independent of the concentration optimization algorithm in section 4.1.5.

### 4.1.4 Saving meltR.F outputs

meltR.F results can be saved to the disk using the "Save_results" argument.

```{r}
meltR.F(df,
        Save_results = "all")
```

The "file_prefix" argument can be used to add a custom file name to the outputs.

```{r}
meltR.F(df,
        Save_results = "all",
        file_prefix = "Helix_J")
```

This will create three pre-canned outputs. The first output, corresponding to Method 1, is a Van't Hoff plot (Figure 5A). Points represent the K\textsubscript{D} and error from fitting isotherms individually. The red line represents the fit to Equation 6 that provides thermodynamic parameters. The blue line and orange line represents the lower and upper limit of the range of K\textsubscript{D} values included in the fit. The second output, corresponding to Method 2, is a depiction of the global fit, where points represent raw data and red lines represent the global fit (Figure 5B). The third output is a .csv file containing the thermodynamic parameters from each method Figure 5 C. Note that the fit line (red) in Figure 5A does not follow the trend at high temperatures and the helix folding energies from this line are thus unreliable. We will cover refining the fit in section 4.1.5.

![Pre-canned meltR.F outputs.](https://user-images.githubusercontent.com/63312483/165543584-0a15344b-9f9c-4b63-ba60-256a7fe2edd1.svg)

### 4.1.5 Refining meltR.F fits

Two arguments are important for refining meltR.F fits. The first is, "Kd_range", which is the range of K\textsubscript{D}s in nM that will be fit to obtain thermodynamic parameters. By default, the "Kd_range" is set to 10 to 1000. The second is, "Kd_error_quantile", which controls the amount of error that is included in the K\textsubscript{D}s that will be fit to obtain thermodynamic parameters. By default, the "Kd_error_quantile" is 0.25, meaning only the 25% most accurate KDs in the "K_range" will be fit to obtain thermodynamic parameters. 

As a first guess, the K\textsubscript{D} range should start about 10 times less than the fluorophore labeled RNA strand concentration and end at about 10 times more than the fluorophore labeled RNA strand, and the "Kd_error_quantile" should be conservative, near 0.25. After this, the Van't Hoff plot should be inspected, for how well the fit matches the linear range. In general, constraining the "Kd_range" and loostening the "Kd_error_quantile" results in the best fits. For example, a "Kd_range" of 40 to 500 nM and a "Kd_error_quantile" of 0.5 works well for the sample data (Figure 6).

```{r}
meltR.F(df,
        Kd_range = c(40, 500),
        Kd_error_quantile = 0.5,
        Save_results = "all")
```

with an output of:

```{r}
[1] "Van't Hoff"
[1] "accurate Ks = 18"
        Method         H     SE.H         S      SE.S         G       SE.G
1    1 VH plot -60.60235 0.964952 -160.5832  3.057177 -10.79747 0.01834264
2 2 Global fit -59.90460 6.321305 -158.3554 20.049234 -10.79066 0.11245376
   K_error         R   Kd.opt
1 6.182951 0.7797618 6.165553
2 6.182951 0.7797618 6.165553
[1] "Fractional error between methods"
           H          S            G
1 0.01158036 0.01397017 0.0006300081
> 
```

Note two things: (1) the refined fit uses more isotherms than the original iterations of the program, where we were checking the concentration optimization algorithm in Section 4.1.3, with an "accurate Ks" count of 18 instead of 10 to 14. (2) The the thermodynamic parameters are closer to the fits with constrained "low_Ks", with a difference of about 5% in terms of the enthalpy and 1% in terms of the free energy.   

![Refined Van't Hoff plot to determine the dH and dS.](https://user-images.githubusercontent.com/63312483/165550513-87659c1c-0a43-4e79-a9ef-d2b22a53e673.svg)

### 4.1.6 Advanced analysis of meltR.F outputs in R

meltR.F can pass a more extensive output to an object in R. 

```{r}
MeltR.fit = meltR.F(df,
        Kd_range = c(40, 500),
        Kd_error_quantile = 0.5,
        Save_results = "all")
```

The object, "MeltR.fit" is now a list of objects that can be passed to plotting functions.

```{r}
names(MeltR.fit)
```

Which will print the names of the list "MeltR.fit".

```{r}
[1] "VantHoff"                         "K"                               
[3] "VH_method_1_fit"                  "VH_method_2_fit"                 
[5] "Raw_data"                         "First_derivative"                
[7] "Tms"                              "R"                               
[9] "Fractional_error_between_methods"
```

[1] VantHoff: is a data frame containing the duflex formation energies. It can be called using:


```{r}
MeltR.fit$VantHoff
```

[2] K: is a data frame containing the results from fitting each isotherm individually. It can be called using:

```{r}
MeltR.fit$K
```

[3] VH_method_1_fit: is a nls object containing the fit obtained from the Van't Hoff plot. It can be called using:

```{r}
MeltR.fit$VH_method_1_fit
```

[4] VH_method_2_fit: is a nls object containing the fit obtained from the global fit. It can be called using:

```{r}
MeltR.fit$VH_method_2_fit
```

[5] Raw_data: The raw data passed back out of MeltR.F with no modifications. It can be called using:

```{r}
MeltR.fit$Raw_data
```

[6] First derivative: The dirst derivative of each sample. Useful for qualitative comparison of data between conditions. It can be called using:

```{r}
MeltR.fit$First_derivative
```

[7] Tms: The approximate T\textsubscript{m} of each sample obtained from the maximum of the first derivative. Useful for qualitative comparison of data between solution conditions. It can be called using:

```{r}
MeltR.fit$Tms
```

[8] R: the mole ratio of fluorophore and quencher labeled RNA, used in the concentration optimization algorithm. The mole ratio "R" was labeled "X" in the theory section to avoid confusion with the gas constant. It can be called using:

```{r}
MeltR.fit$R
```

[9] Fractional error between Methods: The amount thermodynamic parameters vary between methods. It can be called using:

```{r}
MeltR.fit$Fractional_error_between_methods
```

## 4.2 Fitting absorbance melting curves in MeltR

### 4.2.1 Formatting absorbance data for MeltR

Data should be formatted into a comma separated value (“.csv”) text file with carefully labeled columns (Figure 7). There should be a “Sample” column where numbers or character strings specify what sample the absorbance data was collected on. There should be one buffer blank in the data set for background subtraction. If no blank is available, add data with an absorbance of 0 recorded at each temperature. A “Pathlength” column specifies the pathlength in cm of the cuvette the data was collected in. A “Temperature” column specifies the temperature in <span>&#176;</span>C where the data was recorded. Lastly, a “Absorbance” column specifies the absorbance collected at each temperature. 

![Formatting absorbance melting curve data in a csv in Excel. Temperature is in <span>&#176;</span>C and pathlength is in cm.](https://user-images.githubusercontent.com/63312483/165771699-ce0ea0d0-79b5-40c8-8746-9f67b09be245.PNG)

Two common pitfalls occur when formatting data for MeltR, usually in Excel. The first is incorrect column names, even incorporation of an extra space, so that MeltR cannot recognize data when it is read into R. The second is incorporation of characters into data columns when values are missing. If a data point is missing, leave the cell blank. Do not write something like "NA".

### 4.2.2 Reading data into a R data frame

Comma separated value (“.csv”) data can be read into a R data frame for MeltR using the “read.csv” function that is included in base R.

```{r}
df = read.csv("path/file_name.csv")
```

Data can be checked with the "View" function.

```{r}
View(df)
```

Note, the columns in “df” must be named correctly. If the columns are not named correctly, MeltR cannot recognize them. You can rename the columns of a data frame using:

```{r}
colnames(df) = c("Sample", "Pathlength", "Temperature", "Absorbance")
```

### 4.2.3 Applying meltR.A to obtain thermodynamic parameters

Sample data is included in MeltR and can be loaded into your memory and checked.

```{r}
?df.abs.data
df = df.abs.data
```

In general, it is a good idea to plot your data before you start fitting. First, I want to check the quality of low temperature isotherms. Using the tidyverse.

```{r}
head(df)
ggplot(df, aes(x = Temperature, y = Absorbance, color = factor(Sample))) +
  geom_point() +
  theme_classic()
```

![Absorbance data check.](https://user-images.githubusercontent.com/63312483/165774624-39edee9d-8ab8-4ae7-8852-ace5c222303c.svg)

This data looks good (Figure 8), with a small and consistent blank absorbance and sigmoidal RNA melting curves. If one of the RNA melting curves does not resemble a sigmoid, it should be removed from the data set using:

```{r}
df %>% filter(Sample != "Sample you want to remove")
```

Data can now be fit using meltR.A. First we need to define the nucleic acid sequence using a vector that specifies the RNA or DNA sequences used in the experiment.

```{r}
helix = c("RNA", "CGAAAGGU", "ACCUUUCG")
```

meltR.A requires information from 4 user defined arguments to fit the data: the data frame the absorbance data is stored in, the blank sample, the nucleic acid sequences, and the molecular model. Three molecular models are availible in MeltR, "Monomolecular.2State", "Heteroduplex.2State", and "Homoduplex.2State". The absorbance sample data included in MeltR was collected on a heteroduplex. 

```{r}
?meltR.A
meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State")     
```

The output looks like this.

```{r}
[1] "Individual curves"
  Sample           Ct     H       S     G   Tm
1      2 2.832690e-06 -63.5 -0.1725  -9.9 43.1
2      3 4.718659e-06 -71.2 -0.1969 -10.1 44.6
3      4 7.439046e-06 -64.8 -0.1764 -10.0 46.4
4      5 1.173670e-05 -64.5 -0.1752 -10.1 48.4
5      6 2.049582e-05 -70.9 -0.1952 -10.4 50.0
6      7 3.131297e-05 -74.3 -0.2055 -10.5 51.4
7      8 5.533500e-05 -71.6 -0.1974 -10.4 53.0
8      9 8.990554e-05 -74.8 -0.2069 -10.7 54.8
9     10 1.439218e-04 -76.7 -0.2121 -11.0 57.0
[1] "Summary"
              Method     H SE.H      S SE.S     G SE.G
1  1 individual fits -70.3  4.9 -193.1 14.9 -10.4  0.3
2 2 Tm versus ln[Ct] -60.4  1.2 -163.0  3.8  -9.9  0.1
3       3 Global fit -69.7  0.3 -191.3  1.0 -10.4  0.0
[1] "fractional error between methods"
          H         S          G
1 0.1482036 0.1649616 0.04885993
```

Note that the fractional error between the methods is high, over about 15% in terms of the enthalpy. Section 4.2.4 discusses how to refine the fits by trimming the absorbance baselines. 

### 4.2.4 Saving meltR.A outputs

meltR.F results can be saved to the disk using the "Save_results" argument.

```{r}
meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State",
        Save_results = "all")
```

The "file_prefix" argument can be used to add a custom file name to the outputs

```{r}
meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State",
        Save_results = "all",
        file_prefix = "Helix")
```

This will create seven pre-canned outputs. The first and second outputs are depictions of the fits in method Method 1, with the individual fits (red lines) overlay-ed on raw data (Figure 9A) and normalized data generated by calculating the absorbtivity (extinction coefficient) at each temperature (Figure 9B). The third output is a depiction of the fit in Method 2, with the fit overlay-ed on the T\textsubscript{m} data at each strand concentration (Figure 9C). The fourth and fifth outputs are depictions of the fit in method Method 3, with the global fits (red lines) overlay-ed on raw data (Figure 9D) and normalized data generated by calculating the absorbtivity (extinction coefficient) at each temperature (Figure 9E). The sixth output is a .csv file containing the thermodynamic parameters generated by individual fits (Figure 9F). The seventh output is a .csv file containing the thermodynamic parameters from each method (Figure 9G).

![Pre-canned meltR.A outputs](https://user-images.githubusercontent.com/63312483/165786859-fd0e921c-e943-44cf-ad12-7ffd8dd0c20d.svg)

### 4.2.5 Specifying blanks with meltR.A

In this section, we will describe usage for the blank argument, which controls background subtraction. There is a single blank in the data set we fit in section 4.2.3. The identity of the blank sample specified in the "blank" argument, must match the identity of the blank sample in the data frame. If the blank sample is labeled with a character, the "blank" argument needs to be a character string, "Blank 1" for example. In the example in section 4.2.3, the blank sample is represented by a integer. We can check this using the "str" function in base R.

```{r}
str(df.abs.data)
'data.frame':	1800 obs. of  4 variables:
 $ Sample     : int  1 1 1 1 1 1 1 1 1 1 ...
 $ Pathlength : num  1 1 1 1 1 1 1 1 1 1 ...
 $ Temperature: num  4.78 5.5 6.03 6.52 7.03 7.53 8.02 8.54 9.04 9.54 ...
 $ Absorbance : num  0.004395 -0.00058 0.000778 0.001038 0.001068 ...
 ```

As you can see, the Sample column in the data frame, which is why we specified the blank as:

```{r}
blank = 1
```

You may have already performed the background subtraction on your data. In this case, blanks should be removed from the data set and background subtraction can be turned off in meltR.A by setting the blank argument to "none". There is a data set in MeltR that is already background subtracted.

```{r}
?df.FAM.C.BHQ1.data
meltR.A(data_frame = df.FAM.C.BHQ1.data,
        blank = "none",
        NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
        Mmodel = "Heteroduplex.2State")
```

Thus, setting blank = "none" turns off background subtraction.

If there are multiple blanks in the data set, different blanks can be assigned to different blanks by supplying sample-blank pairs to the "blank" argument using a list of vectors. For example, the df.abs.ROX.data data set, included in MeltR, contains variable concentrations of ROX fluorophore. The ROX absorbance must be subtracted out. We first check the Sample identities.

```{r}
?df.abs.ROX.data
unique(df.abs.ROX.data$Sample)
 [1] "0 uM ROX 0 uM Helix A"     "0 uM ROX 5 uM Helix A"    
 [3] "0.032 uM ROX 5 uM Helix A" "0.1 uM ROX 5 uM Helix A"  
 [5] "0.32 uM ROX 5 uM Helix A"  "3.2 uM ROX 5 uM Helix A"  
 [7] "10 uM ROX 5 uM Helix A"    "32 uM ROX 5 uM Helix A"   
 [9] "0.032 uM ROX 0 uM Helix A" "0.1 uM ROX 0 uM Helix A"  
[11] "0.32 uM ROX 0 uM Helix A"  "3.2 uM ROX 0 uM Helix A"  
[13] "10 uM ROX 0 uM Helix A"    "32 uM ROX 0 uM Helix A"  
```

We then create a list of samples specifying sample-blank pairs. Samples are specified first and blanks are specified second.

```{r}
list.blanks = list(c("0 uM ROX 5 uM Helix A", "0 uM ROX 0 uM Helix A"),
              c("0.032 uM ROX 5 uM Helix A", "0.032 uM ROX 0 uM Helix A"),
              c("0.1 uM ROX 5 uM Helix A", "0.1 uM ROX 0 uM Helix A"),
              c("0.32 uM ROX 5 uM Helix A", "0.32 uM ROX 0 uM Helix A"),
              c("3.2 uM ROX 5 uM Helix A", "3.2 uM ROX 0 uM Helix A"),
              c("10 uM ROX 5 uM Helix A", "10 uM ROX 0 uM Helix A"),
              c("32 uM ROX 5 uM Helix A", "32 uM ROX 0 uM Helix A"))
```

We can then supply the list "list.blanks" to the "blank" argument.

```{r}
meltR.A(df.abs.ROX.data,
        blank = list.blanks,
        NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
        methods = c(TRUE, FALSE, FALSE),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(25, 75))
```

Note: I turned of methods 2 and 3 because the data set is not set up for these methods and it would likely cause the algorithm to fail because the data are inappropriate. Also, baselines are trimmed using the fitTs argument. This baseline trimming is discussed in more detail in section 4.2.7.

### 4.2.6 Specifying extinction coefficients

This section will describe how to specify extinction coefficients using the "NucAcid"" argument. In section 4.2.3, we supplied meltR.A with a nucleic acid sequence, and MeltR used extinction coefficient data to calculate an extinction coefficient at 260 nm, which meltR.A then uses to calculate C\textsubscript{t}. MeltR contains data to calculate extinction coefficients at 300, 295, 290, 285, 280, 275, 270, 265, 260, 255, 250, 245, 240, 235, and 230 nm for RNA. MeltR contains data to calculate extinction coefficients at 260 nm for DNA. Extinction coefficient calculations assume pH 7.0.

You can change the wavelength that extinction coefficients are calculated at using the "wavelength" argument. We can test on the df.abs.ROX.data, even though this data was collected at 260 nm.

```{r}
meltR.A(df.abs.ROX.data,
        wavelength = 280,
        blank = blanks,
        NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
        methods = c(TRUE, FALSE, FALSE),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(25, 75))
```

You can check what values meltR.A is calculating using "calc.extcoeff".

```{r}
calc.extcoeff(c("RNA", "CGAAAGGU", "ACCUUUCG"))
$Total
[1] 115200

$CGAAAGGU
[1] 66120

$ACCUUUCG
[1] 49080
```

Custom molar extinction coefficients can be supplied to the NucAcid condition, if data is not included to calculate extinction coefficients for your specific nucleic acid. This is important for modified nucleic acids or extreme solution conditions. For example, the df.FAM.C.BHQ1.data set was collected on a fluorophore and quencher labeled helix. Thus, I need to use extinction coefficients that include the fluorophore and quencher.

```{r}
?df.FAM.C.BHQ1.data
FAM.BHQ1 = c("Custom", 105860, 82300)
meltR.A(data_frame = df.FAM.C.BHQ1.data,
        blank = "none",
        NucAcid = FAM.BHQ1,
        Mmodel = "Heteroduplex.2State",
        concT = 75)
```

### 4.2.7 Refining meltR.A fits by trimming absorbance baselines.

MeltR approximates absorbance baselines using a line (y = mx+b). Real absorbance data deviates from this linear behavior over wide temperature ranges. MeltR fits are thus improved by trimming the baselines so that baselines better approximate a line. First, the raw data should be inspected for linear baselines one sample at a time using the tidyverse. Then, data are fit with meltR.A using a manually optimized baseline trim. After this, the manually optimized fit can be passed to an automated baseline trimmer, which tests combinations of randomly trimmed baselines and identifies an optimum set of trimmed baselines.

```{r}
ggplot(df %>% filter(Sample == 10), #Only plot sample 10
aes(x = Temperature, y = Absorbance, color = factor(Sample))) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = c(15, 70)) #will add horizontal lines to the plot at 15 and 70
```

Thus, for Sample 2, baselines are approximately linear above 80 <span>&#176;</span>C and below 80 <span>&#176;</span>C (Figure 10).

![Baseline trimming on Sample 10. Vertical lines represent 25 and 80 <span>&#176;</span>C.](https://user-images.githubusercontent.com/63312483/165791992-20b92819-5904-4287-9183-c4f489a1d1e7.svg)

This should be repeated for each sample that contains RNA. After identifying ranges where baselines are linear, ranges can be supplied to the "fitTs" argument as a list of vectors (in the order of the samples in the data set). Note, do not provide a temperature range for the blank. 

```{r}
list.T.range = list(c(15, 70), #Sample 2
                    c(15, 70), #Sample 3
                    c(15, 80), #Sample 4
                    c(15, 85), #Sample 5
                    c(20, 78), #Sample 6
                    c(20, 85), #Sample 7
                    c(20, 85), #Sample 8
                    c(25, 85), #Sample 9
                    c(25, 85)) #Sample 10

fit = meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State",
        Save_results = "all",
        file_prefix = "Trimmed",
        fitTs = list.T.range)
```

Which will produce:

```{r}
[1] "Individual curves"
  Sample           Ct     H       S     G   Tm
1      2 2.913343e-06 -67.2 -0.1850  -9.9 42.4
2      3 4.853010e-06 -71.6 -0.1985 -10.1 44.5
3      4 7.650852e-06 -65.3 -0.1784 -10.0 46.1
4      5 1.207087e-05 -61.9 -0.1676  -9.9 47.9
5      6 2.107938e-05 -66.3 -0.1814 -10.0 49.3
6      7 3.220452e-05 -71.2 -0.1964 -10.3 50.8
7      8 5.691051e-05 -68.9 -0.1894 -10.1 52.3
8      9 9.246535e-05 -68.4 -0.1879 -10.1 53.9
9     10 1.480195e-04 -71.0 -0.1953 -10.4 56.1
[1] "Summary"
              Method     H SE.H      S SE.S     G SE.G Tm_at_0.1mM SE.Tm_at_0.1mM
1  1 individual fits -68.0  3.2 -186.7  9.9 -10.1  0.2    54.15446      3.8393625
2 2 Tm versus ln[Ct] -62.1  2.8 -168.9  8.7  -9.8  0.1    53.76500      3.6805415
3       3 Global fit -67.4  0.3 -184.8  0.9 -10.1  0.0    54.26074      0.3579972
[1] "fractional error between methods"
           H          S    G
1 0.08962025 0.09881569 0.03
[1] "dH and dG are in kcal/mol and dS is in cal/mol/K. Tms are in deg Celsius"
```

Fitting manually trimmed baselines reduces the fractional error between the methods from 15% to 9% in terms of the enthalpy.

### 4.2.8 Auto-trimming baselines using the BLtrimmer 

Manually trimming baselines is time consuming, idiosyncratic, and can result in biases. MeltR provides a function called the baseline trimmer that automates baseline trimming by analyzing an ensemble of baseline combinations For more information, please read section 3.2.6. This section will focus on running the BLTrimmer.

The first thing the BLTrimmer needs is a preliminary meltR.A fit object used to identify the core of the melt region that is not trimmed.

```{r}
fit = meltR.A(data_frame = df.abs.data,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State")
```

We can then run this fit through the BLTrimmer. This will take a minute. 

```{r}
?BLTrimmer
BLTrimmer(fit, Save_results = "all")
```

The raw output will yield something like this:

```{r}
[1] "You are trying to test 1000 baseline combinations"
[1] "Do you think this is possible?"
[1] "Fitting 1000 combinations of 5 different baselines per sample"
  |===========================================================================| 100%[1]
[1] "Using autotrimmed baselines in meltR.A"
  |===========================================================================| 100%[1]
  "Ensemble energies"
              Method     dH          CI95.dH      dS            CI95.dS     dG
1  1 individual fits -68.60 -70.02 to -67.32 -188.21 -192.62 to -184.28 -10.23
2 2 Tm versus ln[Ct] -67.99 -69.96 to -66.05 -186.32 -192.43 to -180.33 -10.20
3       3 Global fit -67.93 -70.11 to -66.03 -186.08 -192.85 to -180.09 -10.22
           CI95.dG  Tm_at_0.1mM CI95.Tm_at_0.1mM
1 -10.29 to -10.16        54.68   54.43 to 54.96
2   -10.3 to -10.1        54.69   54.51 to 54.93
3 -10.32 to -10.12        54.79   54.55 to 55.14
[1] "Fractional error between methods"
  dH  dS  dG
1  1 1.1 0.3
[1] "dH and dG are in kcal/mol and dS is in cal/mol/K. Tms are in deg Celsius"
```

This produces thermodynamic parameters from analyzing an ensemble of 250 optimum trimmed baseline combination. The results are saved as a csv and a PDF. This ensemble is shown in Figure 11. Note that the distribution of the "Method 1 dH agrees with Method 2 dH" plot peaks below 5% indicating that the sequence approximates a two-state melting transition.

![BLTrimmer outputs for assessing the quality of the results are written as a PDF file. Blue lines and points represent the 25% of baseline combinations that produce parameters that exhibit the best consistency across the data set. Each point comes from fitting a combination of baselines with Method 1 or 2. Results from 100 baseline combinations are plotted.](https://user-images.githubusercontent.com/63312483/201776066-03f6889b-403b-47d3-bf8b-e6bd9a52bd4e.svg)

We can test the BLTrimmer on a non-two-state folding sequence and see what happens to the distribution.

```{r}
?df.abs.non2S
fit = meltR.A(data_frame = df.abs.non2S,
              NucAcid = c("RNA", "CGCGUUAUAU", "AUAUUUCGCG"),
              Mmodel = "Heteroduplex.2State")
BLTrimmer(fit, Save_results = "all", file_prefix = "Non-two_state")
```

Now the distribution for the "Method 1 dH agrees with Method 2 dH" plot peaks above 15% indicating that this is a non-two-state melting transition (Figure 12).

![BLTrimmer analysis of a non-two-state melting transition. Blue lines and points represent the 25% of baseline combinations that produce parameters that exhibit the best consistency across the data set. Each point comes from fitting a combination of baselines with Method 1 or 2. Results from 100 baseline combinations are plotted.](https://user-images.githubusercontent.com/63312483/201776071-a4d0e88b-9021-4795-a446-4b74bcfaec46.svg)

The BLTrimmer has many parameters for customizing your own baseline trimming protocol. However, if you change the recommended default parameters, you should check to make sure that baselines are not being trimmed to cause the appearence of two-state folding. The data set provided df.abs.data is provided as a positive control. 

### 4.2.8 Advanced plotting meltR.A outputs using the "tidyverse"

meltR.A can pass a more extensive output to an object in R. This also applies to the BLTrimmer. 

```{r}
MeltR.fit = meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State",
        Save_results = "all",
        file_prefix = "Trimmed",
        fitTs = list.T.range)
```

The object, "MeltR.fit" is now a list of objects that can be passed to plotting functions.

```{r}
names(MeltR.fit)
```

Which will print a list of fit names:

```{r}
 [1] "Summary"           "Method.1.indvfits" "Range"            
 [4] "Derivatives.data"  "Method.1.data"     "Method.1.fit"     
 [7] "Method.2.data"     "Method.2.fit"      "Method.3.data"    
[10] "Method.3.fit"   
```

[1] Summary: A data frame containing the thermodynamic parameters from each method. It can be called using:

```{r}
MeltR.fit$Summary
```

[2] Method.1.indfits: A data frame containint the thermodynamic parameters from the individual fits.It can be called using:

```{r}
MeltR.fit$Method.1.indvfits
```

[3] Fractional error between Method 1, 2, and 3. It can be called using:

```{r}
MeltR.fit$Range
```

[4] Derivatives.data: Contains the first and second derivatives for each sample containing RNA. It can be called and plotted using:

```{r}
MeltR.fit$Derivatives.data
ggplot(MeltR.fit$Derivatives.data,
  aes(x = Temperature,
      y = dA.dT/(Pathlength*Ct),
      color = factor(Sample)))+
  geom_point() +
  theme_classic()
```

[5] Method.1.data: Contains the raw data from method 1 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.1.data
ggplot(MeltR.fit$Method.1.data,
    aes(x = Temperature,
        color = factor(Sample),
        group = factor(Sample))) +
  geom_point(mapping = aes(y = Absorbance/(Pathlength*Ct))) +
  geom_line(mapping = aes(y = Model/(Pathlength*Ct)), color = "black") +
  theme_classic()
```
[6] MeltR.fit$Method.1.fit: A list of nls objects containing the fits obtained from fitting melting curves individually. It can be called using:

```{r}
MeltR.fit$Method.1.fit
```

[7] MeltR.fit$Method.2.data: Contains the raw data from method 1 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.2.data
ggplot(MeltR.fit$Method.2.data, aes(x = lnCt)) +
  geom_point(mapping = aes(y = invT)) +
  geom_line(mapping = aes(y = Model)) +
  theme_classic()
```
[8] Method.2.fit: A nls object containing the fit obtained from fitting the relationship of Tm and Ct. It can be called using:

```{r}
MeltR.fit$Method.2.fit
```

[9] Method 3 data: Contains the raw data from method 2 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.3.data
ggplot(MeltR.fit$Method.3.data,
  aes(x = Temperature, 
      color = factor(Sample), 
      group = factor(Sample)))+
  geom_point(mapping = aes(y = Absorbance/(Pathlength*Ct))) +
  geom_line(mapping = aes(y = Model/(Pathlength*Ct)), color = "black") +
  theme_classic()
```
[10]  Method.3.fit: A nls object containing the fit obtained from fitting the raw data. It can be called using:

```{r}
MeltR.fit$Method.3.fit
```

# 5 References

[1] Biochemistry 1996, 35 (45), 14077–14089.
[2] Methods in Enzymology, 1989; Vol. 180, pp 304-325.
[3] Biochemistry 1998, 37 (42), 14719–14735.
[4] Welcome | R for Data Science https://r4ds.had.co.nz/ (accessed 2022 -05 -02).
[5] ggrepel: An R package https://ggrepel.slowkow.com/index.html (accessed 2022 -05 -02).

