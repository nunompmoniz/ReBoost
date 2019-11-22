# ReBoost
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

An R Package with Boosting and SMOTEBoost implementations for Regression Tasks

Data pre-processing methods such as resampling strategies are the most common approach to tackling imbalanced domain learning tasks.
This package encompasses multiple resampling-based boosting strategies for the task of extreme value prediction/imbalanced regression.

**References**

- Nuno Moniz, Rita Ribeiro, Vitor Cerqueira, Nitesh Chawla (2018). "SMOTEBoost for Regression: Improving the Prediction of Extreme Values", Proceedings of the 2018 IEEE 5th International Conference on Data Science and Advanced Analytics.

**To install from github use the following command lines in R:**

    library(devtools)  # You need to install this package!
    install_github("nunompmoniz/ReBoost",ref="master")

After installation the package can be used loaded by doing:

     library(ReBoost)
