# Automated polymer phase identification
Code for Universal Phase Identification of Block Copolymers from Physics-informed Machine Learning 

Xinyi Fang, Elizabeth A. Murphy, Phillip A. Kohl, Youli Li, Craig J. Hawker, Christopher M. Bates, and Mengyang Gu (2024)

This software is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 2, April 2013.

This repository contains the data and code necessary to reproduce the numerical results presented in our research paper. The structure of the repository is organized as follows:

`data/`: This directory holds the copolymer information and SAXS data for all diblock copolymer samples analyzed in the study.

`src/`: Contains the necessary C++ functions for high-performance computing. These functions are interfaced with R using Rcpp and RcppEigen.

`functions/`: Contains R functions written for data analysis.

`diblock_hetero_on_new_chem.R`: Script for predicting phases of copolymers with new chemical identities.

`diblock_hetero_mixed_train_test.R`: Script for predicting phases of copolymers using a mixed dataset of all four groups with 5-fold cross-validation.

`plots_with_all_samples.R`: Script for generating additional plots.
