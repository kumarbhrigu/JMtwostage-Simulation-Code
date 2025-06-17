# JMtwostage-Simulation-Code

This repository contains the R code used for the simulation study in our manuscript titled:  
**"A Two-stage Joint Modeling Approach to Handle Incomplete Time-Dependent Markers in Survival Data through Inverse Probability Weighting and Multiple Imputation."**

The code implements a two-stage method that combines **multiple imputation (MI)** and **inverse probability weighting (IPW)** to address missing longitudinal biomarker data in joint modeling of longitudinal and survival outcomes. The simulations evaluate the performance of the proposed method under various missing data mechanisms (MCAR, MAR, MNAR) and sample sizes.

## 🔧 Installation

You can install the accompanying R package `JMtwostage` directly from GitHub using the [`remotes`](https://cran.r-project.org/package=remotes) package:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install JMtwostage from GitHub
remotes::install_github("kumarbhrigu/JMtwostage")
```
