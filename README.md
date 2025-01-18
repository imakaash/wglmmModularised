
# wglmmModularised

## Modularised Survey-weighted Generalized Linear Mixed Models Allowing for Crossed Effects and Unit Specific Survey Weights

This R package as part of my master thesis submission at University of Trier for sose24 implements a modularised approach to fitting survey-weighted Generalized Linear Mixed Models (GLMMs) that can handle crossed effects and unit-specific survey weights. It provides a flexible framework for analyzing complex survey data with multilevel structures.

## Key Features:

-   Supports various types of GLMMs
-   Handles crossed random effects
-   Incorporates unit-specific survey weights
-   Modular design for flexibility and extensibility

## Installation

To install this package, you need access to the private GitLab repository. Use the following R code:

r

`remotes::install_gitlab("s4akyada/wglmmModularised")` 

**Note:** Your system must be set up for C++ compilation. For Windows, visit [Rtools](https://cran.r-project.org/bin/windows/Rtools/). For Linux and macOS, consult your specific distribution's documentation.This package is designed to provide researchers and data scientists with a powerful tool for analyzing complex survey data, offering advanced modeling capabilities while accounting for survey design features.
