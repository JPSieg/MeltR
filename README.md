# MeltR

[![DOI](https://zenodo.org/badge/262172873.svg)](https://zenodo.org/badge/latestdoi/262172873)

Automated fitting of RNA/DNA absorbance melting curves and fluorescence binding isotherms in R

Jacob P. Sieg<sup>[1],[2]</sup>, Philip C. Bevilacqua<sup>[1],[2],[3]</sup>

<sup>[1]</sup>Department of Chemistry, Pennsylvania State University, University Park, PA 16802.

<sup>[2]</sup>Center for RNA Molecular Biology, Pennsylvania State University, University Park, PA 16802.

<sup>[3]</sup>Department of Biochemistry and Molecular Biology, Pennsylvania State University, University Park, PA 16802.

Last update, 05/16/2022

# 1 Overview

MeltR is a R package, written by Jacob Sieg, that fits nucleic acid folding data to molecular models to obtain thermodynamic parameters. MeltR automates the trivial but time-consuming tasks associated with fitting nucleic acids thermodenaturation data, leading to facile conversion of raw data into useful thermodynamic parameters.

MeltR was inspired by Meltwin.<sup>[1]</sup> MeltR and Meltwin have the same utility: easy and consistent fitting to obtain thermodynamic parameters. The main drawback of MeltR is that it is ran from your R console, whereas Meltwin has a graphical-user-interface. However, the MeltR syntax is not complicated, and MeltR has other advantages: (1) A current versions of MeltR can be downloaded from GitHub by entering two lines of code in your R console, whereas Meltwin has been out of support for years. (2) MeltR supports fitting fluorescence binding isotherms to obtain thermodynamic parameters. (3) MeltR can be ran in bulk. (4) MeltR can be run on any oporating system that is compatible with R. (5) Anecdotally, MeltR is more robust than Meltwin and requires less input from the user.

The core of MeltR is the “meltR.A” function for fitting absorbance melting curves and the “meltR.F” function for fitting fluorescence binding isotherms.

The creators of MeltR are interested in improving MeltR by adding new molecular models. Areas of interest include multistate models and models to describe G-quadruplex RNA. Please email Jacob Sieg at jus841@psu.edu for suggestions.

# 2 Cite MeltR

If you use MeltR to fit fluorescence binding isotherms please cite:

Sieg, J. P.; McKinley, L. N.; Huot, M. J.; Yennawar, N. H.; Bevilacqua, P. C. The Metabolome Weakens RNA Thermodynamic Stability and Strengthens RNA Chemical Stability. Biochemistry 2022. https://doi.org/10.1021/acs.biochem.2c00488.

If you use MeltR to fit absorbance melting curves please cite:

Sieg, J. P.; Arteaga, S. J.; Znosko, B. M.; Bevilacqua, P. C. MeltR Software Provides Facile Determination of Nucleic Acid Thermodynamics. Biophysical Reports 2023, 3 (2), 100101. https://doi.org/10.1016/j.bpr.2023.100101.

# 3 MeltR installation

MeltR is written in R. You can find instructions for installing R here: "https://www.r-project.org/"

MeltR can be installed from "https://github.com/JPSieg/MeltR" by typing two lines into your R console. 

Note, you may need to run R as an administrator for installing devtools and MeltR on on Windows. Likewise, on Mac, the you may need to install the Xcode toolbox from the app store to install devtools.

```{r}
install.packages("devtools")
devtools::install_github("JPSieg/MeltR")
```

MeltR is written entirely in base R and requires no other dependencies. However, Rstudio, the tidyverse, and ggrepel, are recommended software and packages for running MeltR. You can find instructions for installing Rstudio here: "https://www.rstudio.com/products/rstudio/download/". The free version is excellent. Tidyverse and ggrepel can be installed from your R console.

```{r}
install.packages("tidyverse")
install.packages("ggrepel")
```

# 4 Manual

The manual, which contains theory and step-by-step instructions for running MeltR, can be downloaded at this address:

https://github.com/JPSieg/MeltR/blob/master/Manual.pdf

# 5 Video tutorials

1.) Basic workflow: fitting fluorescence binding isotherms with MeltR. https://youtu.be/VNSClpCylDA

2.) Basic workflow: fitting absorbance melting curves with MeltR: https://youtu.be/KkxBFuwJgbE


