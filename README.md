


# MeltR
Automated fitting of RNA/DNA absorbance melting curves and fluorescence binding isotherms in R

# 1 Overview

MeltR is a R package that fits nucleic acid folding data to molecular models to obtain thermodynamic parameters. MeltR automates the trivial but time-consuming tasks associated with fitting nucleic acids thermodenaturation data, and uses “nls” to fit the data to molecular models, leading to facile conversion of raw data into useful thermodynamic parameters.

MeltR was inspired by Meltwin. MeltR and Meltwin have the same utility: easy and consistent fitting to obtain thermodynamic parameters. The main drawback of MeltR is that it is ran from your R console, whereas Meltwin has a graphical-user-interface. However, the MeltR syntax is not complicated, and MeltR has other advantages: (1) A current versions of MeltR can be downloaded from GitHub by entering two lines of code in your R consol, whereas Meltwin has been out of support for years. (2) MeltR supports fitting fluorescence binding isotherms to obtain thermodynamic parameters. (3) MeltR can be ran in bulk. (5) Anecdotally, MeltR is more robust than Meltwin and requires less input from the user.

The core of MeltR is the “meltR.A” function for fitting absorbance melting curves and the “meltR.F” function for fitting fluorescence binding isotherms.

# 2 MeltR installation

```{r}
install.packages("devtools")
devtools::install_github("JPSieg/MeltR")
```

# 3 Theory

## Fitting fluorescence isotherms

## Fitting abosorbance melting curve

### Absorbance data preprocessing

### Thermodynamic models for duplex formation

MeltR can obtain thermodynamic parameters from absorbance melting curves from heteroduplex, homoduplex (selfcomplementary), and monomolecular self-structured DNA and RNA melting curves. Thermo-dynamic parameters for helix formation are obtained using a Van't Hoff model:

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

Method 1 fits the absorbtivity (extinction coefficient) as a function of time for each sample individually. Base lines are modeled as a first order linear model for the absorbance of the unfolded-single standed and folded-duplex state (Equation 5). 

<img src= "https://render.githubusercontent.com/render/math?math={E = mT %2b b\qquad (5)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}K = mT %2b B \qquad (5)}\qquad (1)}#gh-dark-mode-only">

The absorbtivity of each sample as a function of temperature is a function of the fraction of RNA in the folded-duplex state (DS), as a function of temperature f(T).

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})f(T) + (m_{SS}T %2b b_{SS})(1-f(T))\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})f(T) + (m_{SS}T %2b b_{SS})(1-f(T)) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

f(t) is variable, calculated by the analytic solution of the binding constant. MelR uses Equation 7 for heteroduplexes, Equation 8 for homoduplexes, and Equation 8 monomolecular self-structured RNA. 

Thus, method 1 fits absorbtivity versus temperature is fit to equations 9, 10, and 11 to determine for  heteroduplexes, homoduplexes, and monomolecular self-structured RNA respectively. 

## Fitting fluorescence binding isotherms

# 4 Calculating extinction coefficients in MeltR

# 5 Fitting Absorbance Melting Curves in MeltR

## 5.1 Formatting absorbance data for MeltR

## 5.2 Reading data into an R data frame

## 5.3 Applying “meltR.A” to obtain thermodynamic parameters

## 5.4 Saving “meltR.A” outputs

# 6 Fitting Fluorescence Binding Isotherms in MeltR

## 6.1 Formatting fluorescence data for MeltR

## 6.2 Reading data into R
## 6.3 Applying “meltR.F” to obtain thermodynamic parameters	
## 6.4 Saving “meltR.F” outputs

# Advanced plotting MeltR outputs in ggplot2	10

#Appendix

## Appendix A Utility functions for reformatting OLIS data	10

### Reading a “.o3a” formatted file	10
### Pulling data from a single wavelength from multiple “.o3a” files	10
