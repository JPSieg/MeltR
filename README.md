


# MeltR
Automated fitting of RNA/DNA absorbance melting curves and fluorescence binding isotherms in R

# 1 Overview

MeltR is a R package that fits nucleic acid folding data to molecular models to obtain thermodynamic parameters. MeltR automates the trivial but time-consuming tasks associated with fitting nucleic acids thermodenaturation data, leading to facile conversion of raw data into useful thermodynamic parameters.

MeltR was inspired by Meltwin. MeltR and Meltwin have the same utility: easy and consistent fitting to obtain thermodynamic parameters. The main drawback of MeltR is that it is ran from your R console, whereas Meltwin has a graphical-user-interface. However, the MeltR syntax is not complicated, and MeltR has other advantages: (1) A current versions of MeltR can be downloaded from GitHub by entering two lines of code in your R consol, whereas Meltwin has been out of support for years. (2) MeltR supports fitting fluorescence binding isotherms to obtain thermodynamic parameters. (3) MeltR can be ran in bulk. (4) Anecdotally, MeltR is more robust than Meltwin and requires less input from the user.

The core of MeltR is the “meltR.A” function for fitting absorbance melting curves and the “meltR.F” function for fitting fluorescence binding isotherms.

# 2 MeltR installation

```{r}
install.packages("devtools")
devtools::install_github("JPSieg/MeltR")
```

# 3 Theory

## 3.1 Fitting fluorescence isotherms

MeltR can obtains thermodynamic parameters from fluorescence binding isotherms for heteroduplex DNA and RNA using the non-linear regression function in base R, "nls". In this strategy, a quencher labeled strand is titrated, at about 1 to 1000 nM concentrations, into different wells in a qPCR plate containing a constant concentration of fluorophore labeled strand. The fluorophore labeled strand binds to the quencher labeled strand, reducing the fluorescence emmission (E) and resulting in an apparant fluorescence binding isotherm, where the shape of the curve is determined by Equation 1.

<img src= "https://render.githubusercontent.com/render/math?math={E = F_{max} %2b (F_{min} - F_{max})*F(K_{D}, [A]_{T}, [B]_{T}) \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = F_{max} %2b (F_{min} - F_max)*F(K_{D}, [A]_{T}, [B]_{T}) \qquad (1)}#gh-dark-mode-only">

Where Fmax is the fluorescence emission of unbound fluorophore labeled strand (A), Fmin is the fluorescence emission of the fluorophore labeled strand completely bound to a quencher labeled strand (B), and F(K<sub>D</sub>, [A]<sub>T</sub>, [B]<sub>T</sub>) is the mole fraction of the bound fluorophore labeled strand to the total fluorophore labeled strand. 

<img src= "https://render.githubusercontent.com/render/math?math={ F(K_{D}, [A]_{T}, [B]_{T}) = \frac{[AB]}{[A]_{T}} \qquad (2) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} F(K_{D}, [A]_{T}, [B]_{T}) = \frac{[AB]}{[A]_{T}} \qquad (2) }#gh-dark-mode-only">

F(K<sub>D</sub>, [A]<sub>T</sub>, [B]<sub>T</sub>)  is a function controled experimental variables, the total fluorophore labled strandd ([A]<sub>T</sub>) and the total quenceher labeled strand ([B]<sub>T</sub>), and the K<sub>D</sub> is given by the the expression:

<img src= "https://render.githubusercontent.com/render/math?math={K_{D} = \frac{[A][B]}{[AB]} \qquad (3) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} K_{D} = \frac{[A][B]}{[AB]} \qquad (3)}#gh-dark-mode-only">

F(K<sub>D</sub>, [A]<sub>T</sub>, [B]<sub>T</sub>) can be determined by solving the K<sub>D</sub> expression to obtain:

<img src= "https://render.githubusercontent.com/render/math?math={F(K_{D}, [A]_{T}, [B]_{T})  = \frac{(K_{D}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(K_{D}%2b[A]_{T}%2bB]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}} \qquad (4) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} F(Kd, [A]T, [B]T) = \frac{(K_{D}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(K_{D}%2b[A]_{T}%2bB]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}  \qquad (4)}#gh-dark-mode-only">

Thus, the K<sub>D</sub> at each temperature was determined by fitting isotherms at each temperature to equation 5, obtained by plugging equation 4 into equation 1.

<img src= "https://render.githubusercontent.com/render/math?math={E = F_{max} %2b (F_{min} - F_{max})\frac{(K_{D}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(K_{D}%2b[A]_{T}%2bB]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}} \qquad (5) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = F_{max} %2b (F_{min} - F_{max})\frac{(K_{D}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(K_{D}%2b[A]_{T}%2b[B]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}  \qquad (5)}#gh-dark-mode-only">

Thermodynamic parameters for helix formation were extracted by the Van't hoff relationship (Equation 6) between the K<sub>D</sub> the temperature using the two methods described below.

<img src= "https://render.githubusercontent.com/render/math?math={ln(K_{D}) = \frac{dS}{R} - \frac{dH}{RT} \qquad (6) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}ln(K_{D}) = \frac{dS}{R} - \frac{dH}{RT} \qquad (6)}#gh-dark-mode-only">

For a single experiment, all samples should be diluted from the same fluorophore and quencher stocks, so that MeltR can deal with systematic uncertainty in strand concentrations. 

### 3.1.1 Fluorescence data preprocessing

MeltR performs the following data preprocessing steps before fitting:

1.) The fluorophore labeled strand concentration is optomized using the concentration optimization algorithm, described in more detail in section 3.1.3.

2.) The Tm of each sample is estimated using the first derivative of fluorescence emission as a function of temperature. First derivatives are calculated using polynomial regression. The data are fit too a 20th order polynomial to aproximate the data. Then, the first derivative of the polynomial are recorded using calculus. The approximate Tm (where 50% of the nucleic acid is single stranded) is calculated by finding the maximum of the first derivative curve to a precision of less than 0.1 degC. These Tms are unreliable because fluorecence baselines vary wildly with temperature. However, Tms are useful for qualitative comparison of stability between conditions.  

2.) Isotherms are fit to equation x to determine Kd and error in the Kd at each temperature using nls in base R. Initial values for Fmax, Fmin, and Kd are provided by the user. 

3.) Kds are filtered by Kd magnitude and error according to the users specifications to determine which isotherms are most reliable. First, Kds outside of a user specified range (1 to 250 nM by default) are thrown out. Second, Kds are ranked by the error in the Kd. Kds that are below a user specified error quantile are thrown out. The default Kd error quantile file is 0.25, meaning the algorithem with keep the to 25% most accurate Kds.  

4.) The most reliable Kds and temperatures are passed to Method 1 to make Van't Hoff plots.

5.) Fluorescence data from the most reliable Kds and temperatures are passed to Method 2 for global fitting. 

### Concentration optimization algorithm.

We have determined that the accuracy of fit results is highly dependent on errors in the concentration of fluorophore and quencher RNA strands in stock solutions. Fit accuracy is dependent on the mole ratio of fluorophore and quencher labeled RNA in stock solutions, but not dependent on the total magnitude of both fluorophore and quencher labeled RNA concentrations in the stocks. Thus, the concentration optimization algorithm in MeltR does not need to find exact concentrations, just the mole ratio (R) of fluorophore and quencher labeled RNA in samples that should have equal concentrations of fluorophore and quencher, for example 200 nM FAM-RNA and 200 nM RNA-BHQ1.

MeltR uses a fluorecence binding isotherm from a temperature where the Kd is more than 10 times less than the fluorophore labeled strand concentration, where binding to the quencher strand is over-determined. Under these conditions, the shape of the curve will be independent of the Kd, and the isotherm can be used as a Job plot. For example, at a 200 nM fluorophore labeled strand concentration, the shape of the binding curve will be indendent of Kd if the Kd is less than 10 nM. The curve will resemble a hockey stick (Figure 1A) composed of two straint lines. The first line, where [A]T > [B]T, will decrease as the [B]T is increased and the fluorophore labeled strand is saturated. The second line, where [A]T < [B]T, will be constant as the [B]T is increased because the fluorphore labeled strand is saturated. The intersection of the first and second line will occure at:

![Job_plot](https://user-images.githubusercontent.com/63312483/165352390-3e58a2df-1920-4cb4-9805-dee99d67eb6a.svg)

### Figure 1 Concentration optimization algorithm

<img src= "https://render.githubusercontent.com/render/math?math={ [A]_{T} = [B]_{T}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} [A]_{T} = [B]_{T} \qquad (1)}#gh-dark-mode-only">

The absolute concentration of [A]T and [B]T cannot be known with precicion. However, the mole ratio of the error (R) in the fluorophore and quencher stocks can be estimated from this data point at equation x.

<img src= "https://render.githubusercontent.com/render/math?math={ \frac{[A]_{T-estimated}}{R} = [B]_{T}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  \frac{[A]_{T-estimated}}{R} = [B]_{T} \qquad (1)}#gh-dark-mode-only">

Where [A]T-estimated is the estimated total A concentration, or what it should be based on the concentration of the stock, and R is given by equation x.

<img src= "https://render.githubusercontent.com/render/math?math={ R = \frac{[B]_{T}}{[A]_{T-estimated}}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  R = \frac{[B]_{T}}{[A]_{T-estimated}} \qquad (1)}#gh-dark-mode-only">

MeltR fits an overdetermined isotherm to a binding curve (selected by the user but by default the isotherm collected at the lowest temperature) to a modified version of Equation x to determine R (Figure 1 B).

<img src= "https://render.githubusercontent.com/render/math?math={E = F_{max} %2b (F_{min} - F_{max})*\frac{(K_{D}%2b[A]_{T-estimated}%2b[B]_{T}R) - \sqrt{{(K_{D}%2b[A]_{T-estimated}%2bB]_{T}R)}^2 - 4[A]_{T-estimated}[B]_{T}R}}{2[A]_{T-estimated}}  }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = F_{max} %2b (F_{min} - F_{max})*\frac{(K_{D}%2b[A]_{T-estimated}%2b[B]_{T}R) - \sqrt{{(K_{D}%2b[A]_{T-estimated}%2bB]_{T}R)}^2 - 4[A]_{T-estimated}[B]_{T}R}}{2[A]_{T-estimated}}  \qquad (1)}#gh-dark-mode-only">

By default, the fit to determine R allows R and KD to float. The user should also use an argument called "low_K", to set the Kd in the optimization fit to several KDs that are more than 10 times less than the fluorophore labeled strand concentrations. They should then inspect the R from several iterations of the optimization algorithm set to different values to make sure it is similar to the iteration that allows the Kd to float.

[A]T is then corrected with Equation x.

<img src= "https://render.githubusercontent.com/render/math?math={ [A]_{T}  = \frac{[A]_{T-estimated}}{R}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} [A]_{T} = \frac{[A]_{T-estimated}}{R} \qquad (1)}#gh-dark-mode-only">

### Method 1 Van't Hoff plot

Method 1 fits Kds that were passed from the preprocessing steps to equation x. The free energy of helix formation at 37 degC (dG) was calculated from the dH and dS values provided by the fit.

The standard error of the dG was propagated from the standartd error of dH and the dS from the fit, and the covariation between the two values.

<img src= "https://render.githubusercontent.com/render/math?math={ SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dS}}{dS})}^2 %2b 2*310.15*\frac{Covar_{dH, dS}}{dH*dS}} \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dS}}{dS})}^2 %2b 2*310.15*\frac{Covar_{dH, dS}}{dH*dS}}  \qquad (1) }#gh-dark-mode-only">

Helix formation energies are traditionally reported in terms of the association constant.

<img src= "https://render.githubusercontent.com/render/math?math={ K = \frac{[AB]}{[A][B] = 1/K_{D}} \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} K = \frac{[AB]}{[A][B] = 1/K_{D}} \qquad (1)}#gh-dark-mode-only">

Thus, dHs and dSs reported by MeltR are also reported in terms of the association constant, which is obtained by multiplying helix formation energies from fitting Kds by negative one. 

### Method 2 Global fit

Method 2 globally fits fluorescence data that were passed from the preprocessing to equation x, obtained by plugging equation y into equation z. Fmax and Fmin was allowed to float between temperatures and dH and dS were constant.

<img src= "https://render.githubusercontent.com/render/math?math={E = F_{max} %2b (F_{min} - F_{max})*\frac{(\exp{\frac{dS}{R} - \frac{dH}{RT}}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(\exp{\frac{dS}{R} - \frac{dH}{RT}}%2b[A]_{T}%2bB]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}  }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = F_{max} %2b (F_{min} - F_{max})*\frac{(\exp{\frac{dS}{R} - \frac{dH}{RT}}%2b[A]_{T}%2b[B]_{T}) - \sqrt{{(\exp{\frac{dS}{R} - \frac{dH}{RT}}%2b[A]_{T}%2bB]_{T})}^2 - 4[A]_{T}[B]_{T}}}{2[A]_{T}}  \qquad (1)}#gh-dark-mode-only">

Helix folding energies 

## Fitting abosorbance melting curves

### Absorbance data preprocessing

MeltR can obtain thermodynamic parameters from absorbance melting curves for heteroduplex, homoduplex (selfcomplementary), and monomolecular self-structured DNA and RNA using the non-linear regression function in base R, "nls".

MeltR performs the data preprocessing steps before fitting:

1.) RNA and DNA extinction coefficients are calculated from the sequence using the data from (Methods in Enzymology, 1989; Vol. 180, pp 304-325).

2.) The background absorbance of a user specified blank is scaled to the pathlength of the cuvette and subtracted from each curve.

3.) The total strand concentration (Ct) is calculated at a user specified temperature (Default = 90 degC).

4.) High and low temperature data are trimmed (to the users specification) to ensure linear baselines.

5.) First and second derivatives are taken using polynomial regression. First, the data are fit too a 20th order polynomial to aproximate the data. Then, the first and second derivatives of the polynomial are recorded using calculus. The approximate T0.5 (the approximate melting temperature Tm where 50% of the nucleic acid is single stranded) is calculated by finding the maximum of the first derivative curve and the T0.75 (the approximate temperature where 75% of the nucleic acid is single stranded) is calculated by finding the minimum of the first derivative (Figure X), to a precision of less than 0.1 degC.

![alt text](https://user-images.githubusercontent.com/63312483/164780663-09a53ad4-0699-4e29-a31a-c79eaa1b2669.svg "Employee Data title")

6.) Initial parameter estimates are calculated for each curves. The initial values for slopes and intercepts of the baselines are estimated by fitting absorbance values that are greater than the 75th quantile for the uper baselines and fitting aborbance values theat are lower than the 25th quantile for the lower baseline, to y = mx + b. Initial values for the enthalpy are determined using the T0.5 and T0.75 (in Kelvin) from first and second derivative curves. MeltR uses equation x, equation y, and equation zm, for heteroduplex, homoduplex, and monomolecular self-structured DNA and RNA melting curves.   

<img src= "https://render.githubusercontent.com/render/math?math={dH = -0.007*(\frac{1}{T0.5} - \frac{1}{T0.75}) \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}dH = -0.007*(\frac{1}{T0.5} - \frac{1}{T0.75})\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={dH = -0.0044*(\frac{1}{T0.5} - \frac{1}{T0.75}) \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}dH = -0.0044*(\frac{1}{T0.5} - \frac{1}{T0.75})\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={dH = -0.0032*(\frac{1}{T0.5} - \frac{1}{T0.75}) \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}dH = -0.0032*(\frac{1}{T0.5} - \frac{1}{T0.75})\qquad (1)}#gh-dark-mode-only">

The initial Tm is estimated from the first derivative of the melting curve, determined in step 5.

### Thermodynamic models for duplex formation

Thermo-dynamic parameters for helix formation are obtained using a Van't Hoff model:

<img src= "https://render.githubusercontent.com/render/math?math={lnK = \frac{dS}{R} - \frac{dH}{RT}\qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}lnK = \frac{dS}{R} - \frac{dH}{RT}\qquad (1)}#gh-dark-mode-only">

where dS is the entropy change, dH is the enthalpy change, R is the gas constant in kcal/mol, T is the temperature in Kelvin, and K is the equillibrium constant given by equation 2, equation 3, and equation 4 for heteroduplexes, homoduplexes, and monomolecular self-structured RNA respectively.

<img src= "https://render.githubusercontent.com/render/math?math={K = \frac{[AB]}{[A][B]}\qquad (2)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K = \frac{[AB]}{[A][B]}\qquad (2)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={K = \frac{[AA]}{[A]^2}\qquad (3)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K = \frac{[AB]}{[A]^2}\qquad (3)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={K = \frac{[F]}{[U]^2}\qquad (4)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K = \frac{[F]}{[U]^2}\qquad (4)}\qquad (1)}#gh-dark-mode-only">

For equation 2 and equation 3, [A] and [B] are the concentration of different strands, [AB] is the concentration of strand [A] in a duplex with strand [B], and [AA] is the concentration of a self complementary strand A in a duplex with another self complementary strand A. For equation 4, [F] is the concentration of a monomolecular self-structured RNA in the folded state and [U] is the concentration of a monomolecular sel-structured RNA in the unfolded state.

MeltR uses three methods based on the Van't Hoff equation to calculate thermodynamic parameters: (1) fitting melting curves individually, (2) fitting the thermodenaturation point as a function of temperature, and (3) Global fitting melting curves.

#### Method (1) fitting melting curves individually

Method 1 fits the absorbtion as a function of temperature for each sample individually. Base lines are modeled as a first order linear model for the absorbance of the unfolded-single standed and folded-duplex state (Equation 5). 

<img src= "https://render.githubusercontent.com/render/math?math={E = mT %2b b\qquad (5)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K = mT %2b B \qquad (5)}\qquad (1)}#gh-dark-mode-only">

The absrorbtion of each sample as a function of temperature is a function of the fraction of RNA in the folded-duplex state (DS), as a function of temperature f(T).

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})f(T) + (m_{SS}T %2b b_{SS})(1-f(T))\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})f(T) + (m_{SS}T %2b b_{SS})(1-f(T)) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

f(t) is variable, calculated by the analytic solution of the binding constant. MelR uses Equation 7 for heteroduplexes, Equation 8 for homoduplexes, and Equation 9 monomolecular self-structured RNA.

<img src= "https://render.githubusercontent.com/render/math?math={f(T) = \frac{\frac{2}{K(T)*Ct} %2b 2 - \sqrt{(\frac{2}{K(T)*Ct} %2b 2)^2 - 4}}{2} \qquad (7)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}f(T) = \frac{\frac{2}{K(T)*Ct} %2b 2 - \sqrt{(\frac{2}{K(T)*Ct} %2b 2)^2 - 4}}{2} \qquad (6)}\qquad (7)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={f(T) = \frac{\frac{1}{2*K(T)*Ct} %2b 2 - \sqrt{(\frac{1}{2*K(T)*Ct} %2b 2)^2 - 4}}{2} \qquad (7)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}f(T) = \frac{\frac{2}{K(T)*Ct} %2b 2 - \sqrt{(\frac{2}{K(T)*Ct} %2b 2)^2 - 4}}{2} \qquad (6)}\qquad (7)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={f(T) = \frac{K(T)}{1 %2b K(T)} \qquad (7)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}f(T) = \frac{K(T)}{1 %2b K(T)} \qquad (6)}\qquad (7)}#gh-dark-mode-only">

Where Ct is the total strand concentration. K(T) is the equillibrium constant as a function of temperature, given by Equation 7 for heteroduplexes, Equation 8 for homoduplexes, and Equation 9 monomolecular self-structured RNA.

<img src= "https://render.githubusercontent.com/render/math?math={K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm}) %2b ln(\frac{4}{Ct})}}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))} }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))}}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}#gh-light-mode-only">
<img src= "https://render.githubusercontent.com/render/math?math={\color{white}K(T) = \exp{(\frac{H}{R*Tm} - \frac{1}{Tm})})}#gh-dark-mode-only">

Note, the dS term in K(T) has solved in terms of Tm and Ct, then replaced to increase the ease of estimating initial parameters for non-linear regression and to increase the robustness of the nls algorithm.

Thus, method 1 fits absorbtion versus temperature is fit to equations 9, 10, and 11 to determine thermodynamic prameters for heteroduplexes, homoduplexes, and monomolecular self-structured RNA respectively.

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{4}{Ct}))}*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2*\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{H}{R*Tm} - \frac{1}{Tm} %2b ln(\frac{1}{Ct}))}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2*K(T)*Ct} %2b 2 - \sqrt{(\frac{1}{2*K(T)*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*K(T)*Ct} %2b 2 - \sqrt{(\frac{1}{2*K(T)*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b \exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b \exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b \exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

Free energy at 37 degC (dG) is calculated from the dH and entropy (dS) of helix formation 

<img src= "https://render.githubusercontent.com/render/math?math={ dG = dS - 310.15*dS \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dG = dS - 310.15*dS  \qquad (1) }#gh-dark-mode-only">

The dS of helix formation is calculated from the dH and the Tm.

For heteroduplexes:

<img src= "https://render.githubusercontent.com/render/math?math={ dS = \frac{dH}{Tm} %2b R*ln(\frac{4}{Ct}) \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dS = \frac{dH}{Tm} %2b R*ln(\frac{4}{Ct}) \qquad (1) }#gh-dark-mode-only">

For homoduplexes:

<img src= "https://render.githubusercontent.com/render/math?math={ dS = \frac{dH}{Tm} %2b R*ln(\frac{1}{Ct}) \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dS = \frac{dH}{Tm} %2b R*ln(\frac{1}{Ct}) \qquad (1) }#gh-dark-mode-only">

For monomolecular self-structured RNA/DNA:

<img src= "https://render.githubusercontent.com/render/math?math={ dS = \frac{dH}{Tm} \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dS = \frac{dH}{Tm} \qquad (1) }#gh-dark-mode-only">

Error in the dS and dG is calculated by propagating error in the fit terms dH and dT.

<img src= "https://render.githubusercontent.com/render/math?math={ SE_{dS} = |dS|\sqrt{ {(\frac{SE_{dH}}{dH})}^2 %2b {(\frac{SE_{Tm}}{Tm})}^2 - 2*\frac{Covar_{dH, Tm}}{dH*Tm}} \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  SE_{dS} = |dS|\sqrt{ {(\frac{SE_{dH}}{dH})}^2 %2b {(\frac{SE_{Tm}}{Tm})}^2 - 2*\frac{Covar_{dH, Tm}}{dH*Tm}}  \qquad (1) }#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={ SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dH}}{dH})}^2 %2b {(310.15*\frac{SE_{Tm}}{Tm})}^2 - 2*310.15*\frac{Covar_{dH, Tm}}{dH*Tm}} \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dH}}{dH})}^2 %2b {(310.15*\frac{SE_{Tm}}{Tm})}^2 - 2*310.15*\frac{Covar_{dH, Tm}}{dH*Tm}}  \qquad (1) }#gh-dark-mode-only">

### Method 2 fitting the Tm as a function of Ct

Method 2 fits the relationship between 1/Tm and the total strand concentration Ct. To avoid inaccuracies in Tm determination from first derivative plots or covariation with the H terms in method 1, Tms were determined for Method 2 using a semi-quantitative method. Slopes and intercepts from method 1 used to calculate F(T) at each experimental temperature using the absorbance. 

<img src= "https://render.githubusercontent.com/render/math?math={ F(T) = \frac{A -(m_{SS}T %2b b_{SS})}{(m_{DS}T %2b b_{DS}) %2b (m_{SS}T %2b b_{SS})}\qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}F(T) = \frac{A -(m_{SS}T %2b b_{SS})}{(m_{DS}T %2b b_{DS}) %2b (m_{SS}T %2b b_{SS})}\qquad (1)}#gh-dark-mode-only">

F(T) is approimatelty linear in the range of 0.4 to 0.6. Thus, F(T in {0.4 to 0.6}) was fit with y = mT + b, and solved using y = 0.5 to accurately determine the melting temperature for each Ct. To determine thermodynamic parameters, the relationship between 1/Tm and the total strand concentration was then fit to Equations X, and Z, for heteroduplexes and homoduplexes respectively. The Tm of monomolecular, self-structured RNA is independent of Ct so Method 2 cannot be used.

<img src= "https://render.githubusercontent.com/render/math?math={ \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS - R*log(4)}{dH} \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS - R*log(4)}{dH}   \qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={ \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS}{dH}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}F\frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS}{dH}     \qquad (1)}#gh-dark-mode-only">

Free energy at 37 degC (dG) is calculated from the dH and entropy (dS) of helix formation directly from the fit.

<img src= "https://render.githubusercontent.com/render/math?math={ dG = dS - 310.15*dS \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dG = dS - 310.15*dS  \qquad (1) }#gh-dark-mode-only">

Error in the dG is calculated by propagating error in the fit terms dH and dS.

<img src= "https://render.githubusercontent.com/render/math?math={ SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dS}}{dS})}^2 %2b 2*310.15*\frac{Covar_{dH, dS}}{dH*dS}} \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}  SE_{dG} = \sqrt{ {SE_{dH}}^2 %2b {(310.15*\frac{SE_{dS}}{dS})}^2 %2b 2*310.15*\frac{Covar_{dH, dS}}{dH*dS}}  \qquad (1) }#gh-dark-mode-only">

### Method 3 Global fitting

Method 3 fits all curves to equations x, y, and z, simultaneously in a global fit. In this fit, equations x, y, z, are rearranged to be in terms of the dS instead of the Tm. 

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})})*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b\exp{(\frac{dS}{R} - \frac{dH}{RT})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{dS}{R} - \frac{dH}{RT})}}{1 %2b \exp{(\frac{dS}{R} - \frac{dH}{RT})}})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b\exp{(\frac{dS}{R} - \frac{dH}{RT})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{dS}{R} - \frac{dH}{RT})}}{1 %2b \exp{(\frac{dS}{R} - \frac{dH}{RT})}})\qquad (6)\qquad (1)}#gh-dark-mode-only">

The baselines are allowed to vary but dHs and dTs are constrained to a single value for all curves. For global fitting, the slopes and intercepts of the fits from Method 1 are used as initial parameter estimates for the slopes and intercepts of the global fit, and the average of the dHs and dSs from Method 1 are used as initial parameter estimates for the dH and dS.

The dG and error in the dG is calculated using the same equations as Method 2.

# Running MeltR

In this section, we cover how to use MeltR in your R console. If you have not already, read section X to understand the theory underlying the results of MeltR. This section will cover what to type and how to avoid pitfalls. The most common with MeltR is fitting data that is inconsistent with the underlying model, either a fluorescence isotherm or a absorbance melting curve with a non-standard shape. In this case, data will need to be subsetted prior to the fitting set. While MeltR has no dependencies other than base R 4.1.3, data wrangling and plotting functions in the "tidyverse" package are highly recommended, along with the "ggrepel" package, for data quality checks and filtering. To begin, open a new R session in the proper directory. Load relevant packages into your memory.

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

## Fitting fluorescence binding isotherms

### Formatting fluorescence data for MeltR

Data should be formatted into a comma separated value (“.csv”) text file with carefully labeled columns (Figure 5.1). There should be a “Well” column where numbers or character strings specify what well in a microplate the data came from, a “Reading” column that specifies the reading a data point comes from, a “Temperature” column that specifies the temperature where the data was recorded, a “B” column which specifies an approximate quencher labeled strand concentration in nM, an “A” column which specifies an approximate fluorophore labeled strand concentration in nM, and an “Emission” column containing the fluorescence emission intensity. 

![image](https://user-images.githubusercontent.com/63312483/165355843-40b03fbc-7d89-429e-8252-1b65ea6fced1.png)

Two common pitfalls occur when formatting data for MeltR, usually in Excel. The first is incorrect column names, even incorperation of an extra space, so that MeltR cannot recognize relevant data when it is read into R. The second is incorperation of characters characters into data columns when values are missing. If a data point is missing, leave the cell blank. Do not write something like "NA".

### Reading data into an R data frame

Comma separated value (“.csv”) data can be read into an R data frame for MeltR using the “read.csv” function that is included in base R.

```{r}
df = read.csv("path/file_name.csv")
```

Data can be checked with the "View" function.

```{r}
View(df)
```

Note, the columns in “df” should be named correctly. However, if the columns are not named correctly, MeltR cannot recognize them. You can rename the columns of a data frame using:


```{r}
colnames(df1) = c("Helix", "Well", "Reading", "Temperature", "B", "A", "Emission")
```


### Applying “meltR.F” to obtain thermodynamic parameters

For this tutorial, we will use sample data included in MeltR. First check that it is in your memory and find out what is in the data.

```{r}
?df.fluor.data
df = df.fluor.data
```

In general, it is a good idea to plot your data before you start fitting. First, I want to check the quality of low temperature isotherms. Using the tidyverse and ggrepel.

```{r}
head(df)
  Well Reading Temperature   A    B  Emission Helix  Condition
1   A1       1    23.36537 200 1000 0.1671729     J Monovalent
2   A1       2    23.83693 200 1000 0.1680626     J Monovalent
3   A1       3    24.31056 200 1000 0.1687870     J Monovalent
4   A1       4    24.78221 200 1000 0.1695139     J Monovalent
5   A1       5    25.25150 200 1000 0.1703547     J Monovalent
6   A1       6    25.72722 200 1000 0.1712532     J Monovalent
ggplot(df %>% filter(Reading == 1), aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()
```

The result is a really nice binding isotherm (Figure X). One should check a few more readings by changing the value in the filter command, for example 20, 60, 80, etc... If you observe one or two ovious outliers in the data set, it is reasonable to remove them using the filter.

![Isotherm_example](https://user-images.githubusercontent.com/63312483/165367119-ed4e15df-4c08-46d7-82c6-03b386f61395.svg)


```{r}
df = df %>% filter(!Well %in% c("A1", "C3"))
```

If you are satified with the the data, we can move on to fitting. I will first inspect the help file for meltR.F.

```{r}
?meltR.F
```

The only argument that needs set is data frame.

```{r}
fit = meltR.F(df)
[1] "Van't Hoff"
[1] "accurate Ks = 16"
        Method         H     SE.H         S      SE.S         G       SE.G   K_error         R   Kd.opt
1    1 VH plot -48.37473 1.048500 -121.6157  3.357021 -10.65562 0.01033199 0.2916302 0.7797618 6.165553
2 2 Global fit -48.88546 6.302188 -123.2438 20.157522 -10.66140 0.06564137 0.2916302 0.7797618 6.165553
[1] "Fractional error between methods"
           H        S           G
1 0.01050236 0.013298 0.000542719
```

Note, the concentration optimization algorithm, by default allowing the Kd.opt to float, identified an optimal R of ~0.78. Since the FAM-RNA strand in each well is estimated at 200 nM, the concentration optimization algorith will adjust the concentration to 200/0.78 = 256 nM. We will test the robustness of this estimate by constrianing the Kd.op range of KDs that are more than 10 times less than the FAM-RNA concentraion, using the "low_K" argument. For example, 1, 0.5, 0.1, and 0.05 nM. 

```{r}
> meltR.F(df, low_K = 1)
[1] "Van't Hoff"
[1] "accurate Ks = 14"
        Method         H      SE.H         S      SE.S         G        SE.G  K_error         R Kd.opt
1    1 VH plot -60.00212 0.7271876 -158.3742  2.317695 -10.88235 0.009438972 0.345963 0.7231654      1
2 2 Global fit -60.25007 8.9451153 -159.1598 28.490520 -10.88665 0.120521827 0.345963 0.7231654      1
[1] "Fractional error between methods"
            H           S            G
1 0.004123788 0.004947968 0.0003951875
> meltR.F(df, low_K = 0.5)
[1] "Van't Hoff"
[1] "accurate Ks = 14"
        Method         H      SE.H         S      SE.S         G        SE.G   K_error         R Kd.opt
1    1 VH plot -62.21273 0.6849323 -165.3667  2.179732 -10.92425 0.009800243 0.3552451 0.7153701    0.5
2 2 Global fit -62.39456 9.1846917 -165.9409 29.214906 -10.92800 0.134594527 0.3552451 0.7153701    0.5
[1] "Fractional error between methods"
            H           S            G
1 0.002918411 0.003466106 0.0003430056
> meltR.F(df, low_K = 0.1)
[1] "Van't Hoff"
[1] "accurate Ks = 14"
        Method         H      SE.H         S     SE.S         G        SE.G   K_error         R Kd.opt
1    1 VH plot -63.08832 0.6549918 -168.1096  2.08445 -10.94912 0.009371838 0.3609151 0.7080886    0.1
2 2 Global fit -63.28360 9.3633293 -168.7267 29.78152 -10.95302 0.137643619 0.3609151 0.7080886    0.1
[1] "Fractional error between methods"
            H           S            G
1 0.003090634 0.003664033 0.0003556102
> meltR.F(df, low_K = 0.05)
[1] "Van't Hoff"
[1] "accurate Ks = 14"
        Method         H      SE.H         S     SE.S         G        SE.G   K_error         R Kd.opt
1    1 VH plot -63.21087 0.6507751 -168.4935  2.07103 -10.95260 0.009311504 0.3617203 0.7070892   0.05
2 2 Global fit -63.40791 9.3886281 -169.1162 29.86176 -10.95651 0.138076127 0.3617203 0.7070892   0.05
[1] "Fractional error between methods"
            H           S            G
1 0.003112448 0.003688944 0.0003572172
```

The concentration with constrained Kds optomization algorithm adjusts the FAM-RNA concentration to at most 281, about 10% varience with constrained Kds. Thus, the default concentration optimization is robust.

### Saving meltR.F outputs

meltR.F results can be saved to the disk using the "Save_results" argument.

```{r}
meltR.F(df,
        Save_results = "all")
```

The "file_prefix" argument can be used to add a custom file name to the outputs

```{r}
meltR.F(df,
        Save_results = "all",
        file_prefix = "Helix_J")
```

This will create three pre-canned outputs. The first output, corresponding to Method 1, is a Van't Hoff plot (Figure XA). Points represent the Kd and error from fitting isotherms individually. The red line represents the fit to Equation X that provides thermodynamic parameters. The blue line and orange line represents the lower and upper limit of the range of Kd values included in the fit. The second output, corresponding to Method 2, is a depiction of the global fit, where points represent raw data and red lines represent the global fit (Figure XB). The third output is a .csv file containing the thermodynamic parameters from each method.

![meltR F_outputs](https://user-images.githubusercontent.com/63312483/165543584-0a15344b-9f9c-4b63-ba60-256a7fe2edd1.svg)

### Refining meltR.F fits

Two arguments are important for refining meltR.F fits. The first is, "Kd_range", which is the range of KDs in nM that will be fit to obtain thermodynamic parameters. By default, the "Kd_range" is set to 10 to 1000. The second is, "Kd_error_quantile", which controls the amount of error that is included in the KDs that will be fit to obtain thermodynamic parameters. By default, the "Kd_error_quantile" is 0.25, meaning only the 25% most accurate KDs in the "K_range" will be fit to obtain thermodynamic parameters. 

As a first guess, the Kd range should start about 10 times less than the Fluorophore labeled RNA strand concentration and end at about 10 times more than the Fluorophore labeled RNA strand, and the "Kd_error_quantile" should be conservative, near 0.25. After this, the Van't Hoff plot should be inspected. for how well the fit matches the linear range. The "Kd_range" should be constrained and the "Kd_error_quantile"  should be increased. For example, a "Kd_range" of 40 to 500 nM and a "Kd_error_quantile" of 0.5 works well for the sample data (Figure X).

```{r}
meltR.F(df,
        Kd_range = c(40, 500),
        Kd_error_quantile = 0.5,
        Save_results = "all")
```

![Refined_VH_plot](https://user-images.githubusercontent.com/63312483/165550513-87659c1c-0a43-4e79-a9ef-d2b22a53e673.svg)

### Advanced analysis of meltR.F outputs in R

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
[1] "VantHoff"                         "K"                               
[3] "VH_method_1_fit"                  "VH_method_2_fit"                 
[5] "Raw_data"                         "First_derivative"                
[7] "Tms"                              "R"                               
[9] "Fractional_error_between_methods"
```

[1] VantHoff: is a data frame containing the duflex formation energies. It can be called using:


```{r}
MeltR.fit$VantHoff
        Method         H     SE.H         S      SE.S         G       SE.G  K_error         R   Kd.opt
1    1 VH plot -60.60235 0.964952 -160.5832  3.057177 -10.79747 0.01834264 6.182951 0.7797618 6.165553
2 2 Global fit -59.90460 6.321305 -158.3554 20.049234 -10.79066 0.11245376 6.182951 0.7797618 6.165553
```

[2] K: is a data frame containing the results from fitting each isotherm individually. It can be called using:

```{r}
MeltR.fit$K
    Temperature            K        SE.K     Fmax       Fmin
1      23.36537 126470651.49  43997717.6 1.339861 0.15790009
2      23.83693 123142459.99  42167682.9 1.348603 0.15856410
3      24.31056 119560930.01  41121008.6 1.357213 0.15915429
4      24.78221 116170413.96  39360803.6 1.366129 0.16004967
5      25.25150 112810595.89  38232221.2 1.373145 0.16095413
```

[3] VH_method_1_fit: is a nls object containing the fit obtained from the Van't Hoff plot. It can be called using:

```{r}
MeltR.fit$VH_method_1_fit
Nonlinear regression model
  model: lnK ~ Tmodel(H, S, Temperature)
   data: indvfits.to.fit
       H        S 
-60.6024  -0.1606 
 residual sum-of-squares: 0.04117

Number of iterations to convergence: 1 
Achieved convergence tolerance: 3.877e-07
```

[4] VH_method_2_fit: is a nls object containing the fit obtained from the global fit. It can be called using:

```{r}
MeltR.fit$VH_method_2_fit
Nonlinear regression model
  model: Emission ~ Global(H, S, Fmax, Fmin, Reading, A, B, Temperature)
   data: df.gfit
       H        S    Fmax1    Fmax2    Fmax3    Fmax4    Fmax5    Fmax6    Fmax7 
-59.9046  -0.1584   1.5799   1.5862   1.5919   1.5996   1.6050   1.6132   1.6187 
   Fmax8    Fmax9   Fmax10   Fmax11   Fmax12   Fmax13   Fmax14   Fmax15   Fmax16 
  1.6248   1.6325   1.6384   1.6440   1.6513   1.6586   1.6652   1.6706   1.6782 
  Fmax17   Fmax18    Fmin1    Fmin2    Fmin3    Fmin4    Fmin5    Fmin6    Fmin7 
  1.6847   1.6905   0.2683   0.2697   0.2712   0.2727   0.2749   0.2770   0.2803 
   Fmin8    Fmin9   Fmin10   Fmin11   Fmin12   Fmin13   Fmin14   Fmin15   Fmin16 
  0.2848   0.2887   0.2941   0.3004   0.3068   0.3147   0.3243   0.3386   0.3494 
  Fmin17   Fmin18 
  0.3623   0.3803 
 residual sum-of-squares: 3.212
Number of iterations to convergence: 2 
```

[5] Raw_data: The raw data passed back out of MeltR.F with no modifications. It can be called using:

```{r}
Well Reading Temperature        A    B  Emission Helix  Condition
1     A1       1    23.36537 256.4886 1000 0.1671729     J Monovalent
2     A1       2    23.83693 256.4886 1000 0.1680626     J Monovalent
3     A1       3    24.31056 256.4886 1000 0.1687870     J Monovalent
4     A1       4    24.78221 256.4886 1000 0.1695139     J Monovalent
5     A1       5    25.25150 256.4886 1000 0.1703547     J Monovalent
```

[6] First derivative: The dirst derivative of each sample. Useful for qualitative comparison of data between conditions. It can be called using:

```{r}
MeltR.fit$First_derivative
   Well        A    B       Tm        invT           Ct      lnCt
1    A1 256.4886 1000 48.33123 0.003110601 8.717557e-07 -13.95276
2    A2 256.4886  800 47.22254 0.003121366 6.717557e-07 -14.21337
3    A3 256.4886  600 46.47811 0.003128636 4.717557e-07 -14.56680
4    A4 256.4886  400 45.42207 0.003139007 2.717557e-07 -15.11836
5    A5 256.4886  250 44.47293 0.003148387 1.217557e-07 -15.92125
```

[7] Tms: The approximate Tm of each sample obtained from the maximum of the first derivative. Useful for qualitative comparison of data between solution conditions. It can be called using:

```{r}
MeltR.fit$Tms
   Well        A    B       Tm        invT
1    A1 256.4886 1000 48.33123 0.003110601
2    A2 256.4886  800 47.22254 0.003121366
3    A3 256.4886  600 46.47811 0.003128636
4    A4 256.4886  400 45.42207 0.003139007
5    A5 256.4886  250 44.47293 0.003148387
```

[8] R: the mole ratio of fluorophore and quencher labeled RNA, used in the concentration optimization algorithm. It can be called using:

```{r}
MeltR.fit$R
        R 
0.7797618 
```

[9] Fractional error between Methods: The amount thermodynamic parameters vary between methods. It can be called using:

```{r}
MeltR.fit$Fractional_error_between_methods
           H          S            G
1 0.01158036 0.01397017 0.0006300081
```

## Fitting Absorbance Melting Curves in MeltR

### Formatting absorbance data for MeltR

Data should be formatted into a comma separated value (“.csv”) text file with carefully labeled columns (Figure X). There should be a “Sample” column where numbers or character strings specify what sampe the absorbance data was collected. There should be one buffer blank in the data set for background subtraction. If no blank is availible, add data with an absorbance of 0 recorded at each temperature. A “Pathlength” column specifies the pathlength in cm of the cuvette the data was collected in. A “Temperature” column specifies the temperature where the data was recorded. Lastly, a “Absorbance” column specifies the absorbance collected at each temperature. 

![Figure_4 1](https://user-images.githubusercontent.com/63312483/165771699-ce0ea0d0-79b5-40c8-8746-9f67b09be245.PNG)


Two common pitfalls occur when formatting data for MeltR, usually in Excel. The first is incorrect column names, even incorperation of an extra space, so that MeltR cannot recognize data when it is read into R. The second is incorperation of characters characters into data columns when values are missing. If a data point is missing, leave the cell blank. Do not write something like "NA".

### Reading data into an R data frame


Comma separated value (“.csv”) data can be read into an R data frame for MeltR using the “read.csv” function that is included in base R.

```{r}
df = read.csv("path/file_name.csv")
```

Data can be checked with the "View" function.

```{r}
View(df)
```

Note, the columns in “df” should be named correctly. However, if the columns are not named correctly, MeltR cannot recognize them. You can rename the columns of a data frame using:


```{r}
colnames(df1) = c("Sample", "Pathlength", "Temperature", "Absorbance")
```

### Applying meltR.A to obtain thermodynamic parameters

Sample data is included in MeltR. Data can be loaded into your memory and checked.

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

![Check_absorbance_data](https://user-images.githubusercontent.com/63312483/165774624-39edee9d-8ab8-4ae7-8852-ace5c222303c.svg)

This data looks good, with a small and consistent blank absorbance and sigmoidal RNA melting curves. If one of the RNA melting curves does not resemble a sigmoid, it should be removed from the data set using:

```{r}
df = df %>% filter(Sample != "Sample you want to remove")
```

Data can now be fit using meltR.A. First we need to define the nucleic acid sequence using a vector that specifies the RNA or DNA sequences used in the experiment.

```{r}
helix = c("RNA", "CGAAAGGU", "ACCUUUCG")
```

meltR.A requires information from 4 user defined arguments to fit the data: the data frame the absorbance data is stored in, the blank sample, the nucleic acid sequences, and the molecular model. Three molecular models are availible in MeltR, "Monomolecular.2State", "Heteroduplex.2State", and "Homoduplex.2State".

```{r}
?meltR.A
meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State")
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

Note that the fractional error between the methods is over about 15% in terms of the enthalpy. Section X.X discusses how to refine the fits by trimming the absorbance baselines. 

### Saving meltR.A outputs

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

This will create seven pre-canned outputs. The first and second outputs are depictions of the fits in method Method 1, with the individual fits (red lines) overlayed on raw data (Figure xA) and normalized data generated by calculating the absorbtivity (extinction coefficient) at each temperature (Figure XB). The third output is a depiction of the fit in Method 2, with the fit overlayed on the Tm data at each strand concentration (Figure XC). The fourth and fith outputs are depictions of the fit in method Method 3, with the global fits (red lines) overlayed on raw data (Figure xD) and normalized data generated by calculating the absorbtivity (extinction coefficient) at each temperature (Figure XE). The sixth output is a .csv file containing the thermodynamic parameters generated by individual fits (Figure XF). The seventh output is a .csv file containing the thermodynamic parameters from each method (Figure XG).

![Absorbance_outputs](https://user-images.githubusercontent.com/63312483/165786859-fd0e921c-e943-44cf-ad12-7ffd8dd0c20d.svg)


### Refining MeltR.A fits by trimming fluorescence baselines.

MeltR approximates absorbance baselines using a line (y = mx+b). Real absorbance data deviates from this linear behavior over wide temperature ranges. MeltR fits are thus improved by trimming the baselines so that they better approximate a line. First, the raw data should be inspected for linear baselines one sample at a time using the tidyverse.

```{r}
ggplot(df %>% filter(Sample == 2), aes(x = Temperature, y = Absorbance, color = factor(Sample))) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = c(15, 70))
```

Thus, for Sample 2, baselines are approximately linear above 15 degrees C and below 70 degrees C (Figure X).

![Baseline_trimming](https://user-images.githubusercontent.com/63312483/165791992-20b92819-5904-4287-9183-c4f489a1d1e7.svg)

This should be repeated for each sample that contains RNA). After identifying ranges where baselines are linear, ranges can be supplied to the "fitTs" argument as a list of vectors (in the order of the samples in the data set).

```{r}
list.T.range = list(c(15, 70),
                    c(15, 70),
                    c(15, 80),
                    c(15, 85),
                    c(20, 78),
                    c(20, 85),
                    c(20, 85),
                    c(25, 85),
                    c(25, 85))

meltR.A(data_frame = df,
        blank = 1,
        NucAcid = helix,
        Mmodel = "Heteroduplex.2State",
        Save_results = "all",
        file_prefix = "Trimmed",
        fitTs = list.T.range)
[1] "Individual curves"
  Sample           Ct     H       S     G   Tm
1      2 2.832690e-06 -67.2 -0.1849  -9.9 42.4
2      3 4.718659e-06 -71.6 -0.1984 -10.1 44.5
3      4 7.439046e-06 -65.3 -0.1783 -10.0 46.1
4      5 1.173670e-05 -61.9 -0.1676 -10.0 47.9
5      6 2.049582e-05 -66.3 -0.1813 -10.0 49.3
6      7 3.131297e-05 -71.2 -0.1964 -10.3 50.8
7      8 5.533500e-05 -68.9 -0.1894 -10.1 52.3
8      9 8.990554e-05 -68.4 -0.1879 -10.1 53.9
9     10 1.439218e-04 -71.0 -0.1952 -10.4 56.1
[1] "Summary"
              Method     H SE.H      S SE.S     G SE.G
1  1 individual fits -68.0  3.2 -186.6  9.9 -10.1  0.2
2 2 Tm versus ln[Ct] -62.2  1.4 -168.9  4.5  -9.9  0.1
3       3 Global fit -67.4  0.3 -184.7  0.9 -10.1  0.0
[1] "fractional error between methods"
           H          S          G
1 0.08805668 0.09829693 0.01993355
```
The fitting trimmed baselines reduces the fractional error between the methods from 15% to 9%. 

### Advanced plotting meltR.A outputs using the "tidyverse"

meltR.A can pass a more extensive output to an object in R. 

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
 [1] "Summary"           "Method.1.indvfits" "Range"            
 [4] "Derivatives.data"  "Method.1.data"     "Method.1.fit"     
 [7] "Method.2.data"     "Method.2.fit"      "Method.3.data"    
[10] "Method.3.fit"   
```

[1] Summary: A data frame containing the thermodynamic parameters from each method. It can be called using:

```{r}
MeltR.fit$Summary
              Method     H SE.H      S SE.S     G SE.G
1  1 individual fits -68.0  3.2 -186.6  9.9 -10.1  0.2
2 2 Tm versus ln[Ct] -62.2  1.4 -168.9  4.5  -9.9  0.1
3       3 Global fit -67.4  0.3 -184.7  0.9 -10.1  0.0
```

[2] Method.1.indfits: A data frame containint the thermodynamic parameters from the individual fits.It can be called using:

```{r}
MeltR.fit$Method.1.indvfits
  Sample           Ct     H       S     G   Tm
1      2 2.832690e-06 -67.2 -0.1849  -9.9 42.4
2      3 4.718659e-06 -71.6 -0.1984 -10.1 44.5
3      4 7.439046e-06 -65.3 -0.1783 -10.0 46.1
4      5 1.173670e-05 -61.9 -0.1676 -10.0 47.9
5      6 2.049582e-05 -66.3 -0.1813 -10.0 49.3
6      7 3.131297e-05 -71.2 -0.1964 -10.3 50.8
7      8 5.533500e-05 -68.9 -0.1894 -10.1 52.3
8      9 8.990554e-05 -68.4 -0.1879 -10.1 53.9
9     10 1.439218e-04 -71.0 -0.1952 -10.4 56.1
```

[3] Fractional error between Method 1, 2, and 3. It can be called using:

```{r}
MeltR.fit$Range
           H          S          G
1 0.08805668 0.09829693 0.01993355
```

[4] Derivatives.data: Contains the first and second derivatives for each sample containing RNA. It can be called and plotted using:

```{r}
MeltR.fit$Derivatives.data
    Sample Pathlength Temperature Absorbance           Ct         dA.dT        dA.dT2
201      2          1       15.07  0.2717590 2.832690e-06  2.811276e-04 -1.803953e-04
202      2          1       15.58  0.2713165 2.832690e-06  1.967808e-04 -1.463624e-04
203      2          1       16.09  0.2723999 2.832690e-06  1.328851e-04 -1.061979e-04
204      2          1       16.59  0.2719421 2.832690e-06  8.941285e-05 -6.758778e-05
205      2          1       17.10  0.2721405 2.832690e-06  6.375680e-05 -3.442765e-05

ggplot(MeltR.fit$Derivatives.data, aes(x = Temperature, y = dA.dT/(Pathlength*Ct), color = factor(Sample)))+
  geom_point() +
  theme_classic()
```

![First_derivative_data](https://user-images.githubusercontent.com/63312483/165796260-5d6b9070-7b74-4323-b4d5-fafabbcef2e8.svg)

[5] Method.1.data: Contains the raw data from method 1 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.1.data
    Sample Pathlength Temperature Absorbance           Ct       Ext     Model
201      2          1       15.07  0.2717590 2.832690e-06  95936.75 0.2717559
202      2          1       15.58  0.2713165 2.832690e-06  95780.54 0.2717820
203      2          1       16.09  0.2723999 2.832690e-06  96162.99 0.2718099
204      2          1       16.59  0.2719421 2.832690e-06  96001.39 0.2718393
205      2          1       17.10  0.2721405 2.832690e-06  96071.42 0.2718716
ggplot(MeltR.fit$Method.1.data, aes(x = Temperature, color = factor(Sample), group = factor(Sample)))+
  geom_point(mapping = aes(y = Absorbance/(Pathlength*Ct))) +
  geom_line(mapping = aes(y = Model/(Pathlength*Ct)), color = "black") +
  theme_classic()
```

![Method_1_data](https://user-images.githubusercontent.com/63312483/165797375-d9f7ab04-2749-43ff-962f-6e935cf0ace4.svg)

[6] MeltR.fit$Method.1.fit: A list of nls objects containing the fits obtained from fitting melting curves individually. It can be called using:

```{r}
MeltR.fit$Method.1.fit
```

[7] MeltR.fit$Method.2.data: Contains the raw data from method 1 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.2.data
        lnCt       Tm        invT       Model
1 -12.774284 42.45013 0.003168566 0.003165178
2 -12.263986 44.43509 0.003148762 0.003148887
3 -11.808768 46.09716 0.003132369 0.003134353
4 -11.352790 47.81681 0.003115587 0.003119796
5 -10.795290 49.12745 0.003102916 0.003101997
ggplot(MeltR.fit$Method.2.data, aes(x = lnCt)) +
  geom_point(mapping = aes(y = invT)) +
  geom_line(mapping = aes(y = Model)) +
  theme_classic()
```

![Method_2_data](https://user-images.githubusercontent.com/63312483/165801194-97d0934e-ce95-4de9-8ff0-75a236fd6659.svg)

[8] Method.2.fit: A nls object containing the fit obtained from fitting the relationship of Tm and Ct. It can be called using:

```{r}
MeltR.fit$Method.2.fit
Nonlinear regression model
  model: invT ~ TmModel(H, S, lnCt)
   data: Tm_data
       H        S 
-62.2438  -0.1689 
 residual sum-of-squares: 5.535e-11

Number of iterations to convergence: 4 
Achieved convergence tolerance: 2.206e-08
```

[9] Method 3 data: Contains the raw data from method 2 and the model. It can be called and plotted using:

```{r}
MeltR.fit$Method.3.data
    Sample Pathlength Temperature Absorbance           Ct       Ext     Model
201      1          1       15.07  0.2717590 2.832690e-06  95936.75 0.2711654
202      1          1       15.58  0.2713165 2.832690e-06  95780.54 0.2712364
203      1          1       16.09  0.2723999 2.832690e-06  96162.99 0.2713091
204      1          1       16.59  0.2719421 2.832690e-06  96001.39 0.2713822
205      1          1       17.10  0.2721405 2.832690e-06  96071.42 0.2714587

ggplot(MeltR.fit$Method.3.data, aes(x = Temperature, color = factor(Sample), group = factor(Sample)))+
  geom_point(mapping = aes(y = Absorbance/(Pathlength*Ct))) +
  geom_line(mapping = aes(y = Model/(Pathlength*Ct)), color = "black") +
  theme_classic()
```
[10]  Method.3.fit: A nls object containing the fit obtained from fitting the raw data. It can be called using:

```{r}
MeltR.fit$Method.3.fit
Nonlinear regression model
  model: Absorbance ~ GModel(H, S, mED, bED, mESS, bESS, Sample, Temperature,     Ct)
   data: gfit_data
         H          S       mED1       mED2       mED3       mED4       mED5 
-6.741e+01 -1.847e-01  1.083e-04  8.979e-05  2.023e-04  2.624e-05 -1.496e-04 
      mED6       mED7       mED8       mED9       bED1       bED2       bED3 
-1.104e-04  2.011e-04  2.073e-04  2.026e-04  2.694e-01  4.634e-01  7.198e-01 
      bED4       bED5       bED6       bED7       bED8       bED9      mESS1 
 1.161e+00  1.016e+00  1.570e+00  5.321e-01  8.743e-01  1.436e+00  3.638e-04 
     mESS2      mESS3      mESS4      mESS5      mESS6      mESS7      mESS8 
 1.862e-04  5.775e-04  5.843e-04  4.111e-04  2.897e-04  1.899e-04  1.432e-04 
     mESS9      bESS1      bESS2      bESS3      bESS4      bESS5      bESS6 
 6.605e-05  3.046e-01  5.374e-01  8.315e-01  1.343e+00  1.178e+00  1.829e+00 
     bESS7      bESS8      bESS9 
 6.392e-01  1.053e+00  1.707e+00 
 residual sum-of-squares: 0.001631

Number of iterations to convergence: 2 
Achieved convergence tolerance: 0.0002484
```

# References
