# MAIHDA: Perceived Health Trajectories

This repository contains the R code for the study: "Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA) for Perceived Health Trajectories" (Belachew et al., 2026).

The analysis utilizes Bayesian ordinal multilevel modeling to examine intersectional inequalities in percieved health trajectories from adolescence (age 16) to early adulthood (age 33).

## Project Structure:

  * ### Root Directory
    * `MAIHDA_analysis.R`       - Main analysis script
    * `MAIHDA_data_final.rds`   - Final processed dataset (input)
    * `LICENSE` — MIT License for open-source distribution.
    * `README.md` — Documentation and setup instructions.
  * `outputs/`                - auto generated folder
    * `figures/`            - High-resolution TIFF plots (600 DPI)
    * `diagnostics/`       - MCMC trace plots and LOO-CV results, model validation
    * `Appendix_Full_Ranking.xlsx` - Stratum-specific rankings

## Prerequisites

The script is built using R (version 4.0+). It will automatically check for and install the following libraries:

  * `Modeling:` brms (Bayesian Regression Models using 'Stan')
  * `Tidying/Visualization:` tidyverse, tidybayes, bayesplot, patchwork, ggrepel
  * `Data/Export:` haven, openxlsx, matrixStats

## Hardware Note

The Bayesian models are configured to run on 4 cores in parallel. It is recommended to have at least 8GB of RAM, as the moment_match LOO-CV process is computationally intensive.

## Getting Started

### ***1. Clone the repository:***
         git clone https://github.com/Abate-Belachew/MAIHDA-perceived-health-trajectory.git
### ***2. Place your data:*** 
         Ensure `MAIHDA_data_final.rds` is in the root directory. 
### ***3. Run the analysis:*** 
          Open MAIHDA_analysis.R in RStudio and source the file, or run via terminal:
          source("MAIHDA_analysis.R")

## Key Analysis Steps

  * `Intersectional Strata Construction:` Creates unique IDs combining gender, family type, education, income, and social support.

  * `Bayesian Ordinal Models:` Fits a cumulative logit model with random intercepts for strata and project IDs.

  * `Model Comparison:` Uses Leave-One-Out Cross-Validation (LOO-CV) with moment matching for robust model selection.
  * `Trajectory Analysis:` Identifies the top 5 increasing and decreasing health trajectories across life stages.
  * `Output Generation:` Exports production-ready figures (TIFF, 600 DPI) and full appendix tables.

## Citation

If you use this code or method in your research, please cite:

Belachew, A., et al. (2026). Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA) for Perceived Health Trajectories.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
