


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

## Fitting abosorbance melting curves

### Absorbance data preprocessing

MeltR can obtain thermodynamic parameters from absorbance melting curves from heteroduplex, homoduplex (selfcomplementary), and monomolecular self-structured DNA and RNA melting curves using the non-linear regression function in base R, "nls".

MeltR performs the data preprocessing steps before fitting:

1.) RNA and DNA extinction coefficients are calculated from the sequence using the data from (Methods in Enzymology, 1989; Vol. 180, pp 304-325).

2.) The background absorbance of a user specified blank is scaled to the pathlength of the cuvette and subtracted from each curve.

3.) The total strand concentration (Ct) is calculated at a user specified temperature (Default = 90 degC).

4.) High and low temperature data are trimmed (to the users specification) to ensure linear baselines.

5.) First and second derivatives are taken using polynomial regression. First, the data are fit too a 20th order polynomial to aproximate the data. Then, the first and second derivatives of the polynomial. The approximate T0.5 (the approximate melting temperature Tm where 50% of the nucleic acid is single stranded) is calculated by finding the maximum of the first derivative curve and the T0.75 (the approximate temperature where 75% of the nucleic acid is single stranded) is calculated by finding the minimum of the first derivative (Figure X), to a precision of less than 0.1 degC.


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

<img src= "https://render.githubusercontent.com/render/math?math={ dS = \frac{dH}{Tm} %2b R*ln(\frac{1}{Ct}) \qquad (1) }#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} dS = \frac{dH}{Tm} %2b R*ln(\frac{1}{Ct}) \qquad (1) }#gh-dark-mode-only">



calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(4/Ct)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm)))) }
      calcG = function(H, Tm, Ct){ H - (310.15*((H/(273.15 + Tm)) + (0.0019872*log(4/Ct)))) }
      calcG.SE = function(H, Tm, Ct, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*((H/(273.15 + Tm)) + (0.0019872*log(4/Ct))))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2) }

for homoduplexes:


calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(1/Ct)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))}
      calcG = function(H, Tm, Ct){ H - (310.15*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct)))) }
      calcG.SE = function(H, Tm, Ct, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*((H/(273.15 + Tm)) + (0.0019872*log(1/Ct))))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2)}

For monomolecular self-structured RNA/DNA:

calcS = function(H, Tm){ (H/(273.15 + Tm)) }
      calcS.SE = function(H, Tm, SE.H, SE.Tm, covar){ abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm)))) }
      calcG = function(H, Tm){H - (310*((H/(273.15 + Tm))))}
      calcG.SE = function(H, Tm, SE.H, SE.Tm, covar){ sqrt((SE.H)^2 + (abs(310*(H/(273.15 + Tm)))*(abs(H/(273.13 + Tm))*sqrt((SE.H/H)^2 + (SE.Tm/Tm)^2 - (2*(covar/(H*Tm))))))^2) }

### Method 2 fitting the Tm as a function of Ct

Method 2 fits the relationship between 1/Tm and the total strand concentration Ct. To avoid inaccuracies in Tm determination from first derivative plots or covariation with the H terms in method 1, Tms were determined for Method 2 using a semi-quantitative method. Slopes and intercepts from method 1 used to calculate F(T) at each experimental temperature using the absorbance. 

<img src= "https://render.githubusercontent.com/render/math?math={ F(T) = \frac{A -(m_{SS}T %2b b_{SS})}{(m_{DS}T %2b b_{DS}) %2b (m_{SS}T %2b b_{SS})}\qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}F(T) = \frac{A -(m_{SS}T %2b b_{SS})}{(m_{DS}T %2b b_{DS}) %2b (m_{SS}T %2b b_{SS})}\qquad (1)}#gh-dark-mode-only">

F(T) is approimatelty linear in the range of 0.4 to 0.6. Thus, F(T in {0.4 to 0.6}) was fit with y = mT + b, and solved using y = 0.5 to accurately determine the melting temperature for each Ct. To determine thermodynamic parameters, the relationship between 1/Tm and the total strand concentration was then fit to Equations X, and Z, for heteroduplexes and homoduplexes respectively. The Tm of monomolecular, self-structured RNA is independent of Ct so Method 2 cannot be used.

<img src= "https://render.githubusercontent.com/render/math?math={ \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS - R*log(4)}{dH} \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white} \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS - R*log(4)}{dH}   \qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={ \frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS}{dH}  \qquad (1)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}F\frac{1}{Tm} =  \frac{R}{dH}*lnCt) %2b \frac{dS}{dH}     \qquad (1)}#gh-dark-mode-only">

### Method 3 Global fitting

Method 3 fits all curves to equations x, y, and z, simultaneously in a global fit. In this fit, equations x, y, z, are rearranged to be in terms of the dS instead of the Tm. 

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{2}{\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})})*Ct} %2b 2)^2 - 4}}{2} + (m_{SS}T %2b b_{SS})(1-\frac{\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2 - \sqrt{(\frac{1}{2*\exp{(\frac{dS}{R} - \frac{dH}{RT})}*Ct} %2b 2)^2 - 4}}{2}) \qquad (6)}\qquad (1)}#gh-dark-mode-only">

<img src= "https://render.githubusercontent.com/render/math?math={E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b\exp{(\frac{dS}{R} - \frac{dH}{RT})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{dS}{R} - \frac{dH}{RT})}}{1 %2b \exp{(\frac{dS}{R} - \frac{dH}{RT})}})\qquad (6)}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}E = (m_{DS}T %2b b_{DS})\frac{\exp{(\frac{H}{R*Tm} - \frac{1}{Tm})}}{1 %2b\exp{(\frac{dS}{R} - \frac{dH}{RT})}} + (m_{SS}T %2b b_{SS})(1-\frac{\exp{(\frac{dS}{R} - \frac{dH}{RT})}}{1 %2b \exp{(\frac{dS}{R} - \frac{dH}{RT})}})\qquad (6)\qquad (1)}#gh-dark-mode-only">

The baselines are allowed to vary but dHs and dTs are constrained to a single value for all curves. For global fitting, the slopes and intercepts of the fits from Method 1 are used as initial parameter estimates for the slopes and intercepts of the global fit, and the average of the dHs and dSs from Method 1 are used as initial parameter estimates for the dH and dS.


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
