# Interval Truth Model

This repository contains code and supplementary files for the manuscript:

> Kloft, M.,  Siepe, B.S., Heck, D.W. (2024). The Interval Truth Model: A Consensus Model for Continuous Bounded Interval Responses. <https://doi.org/10.31234/osf.io/dzvw2>.

A BiBTeX entry for LaTeX users is:

```BibTeX
@article{Kloft2024,
  title={The Interval Truth Model: A Consensus Model for Continuous Bounded Interval Responses},
  author={Kloft, Matthias and Siepe, Bj{\"o}rn S. and Heck, Daniel W.},
  url={https://doi.org/10.31234/osf.io/dzvw2},
  year={2024},
  note={PsyArXiv preprint}
}
``` 


**Note**
If you encounter the acronym *ITM* in the code or file names, 
it refers to the Interval Truth Model, which was the name of the model in the first version of the manuscript. 
The name was later changed to Interval Consensus Model (ICM) to better reflect its purpose.

## Folder Structure
### `data/`

This folder contains the data used in the empirical example.


### `plots/`

This folder contains plots used in the manuscript and in the supplementary material.

### `sim_results/`

This folder contains the results of the simulation studies. 

### `fits/`

This folder contains saved model fits from the empirical analyses, including fitted Stan models for the verbal quantifiers dataset and sensitivity analyses.

### `src/`

This folder contains the code to reproduce our results, as well as the supplementary files for our simulation studies. 

- `00_download_big_files_osf.R`: Script to download the fit files and simulation results from the OSF storage.
- `00_functions.R`: Auxiliary functions
- `01a_simluation_study_final_v2.Rmd`: Code for the main simulation study. The corresponding `.R`-file was used to run the simulation on a server.
- `01b_sim_visualizations_v2.Qmd`: Visualizations of the results of the main simulation study. The corresponding `.html` file shows a rendered version of the document.
- `01c_simluation_study_link_functions_final.Rmd`: Code for the link function simulation study. The corresponding `.R`-file was used to run the simulation on a server.
- `01d_sim_link_functions_visualizations.Qmd`: Visualizations of the results of the link function simulation study. The corresponding `.html` file shows a rendered version of the document.
- `01e_simulation_development_ilr.Qmd`: Justification and visualization of the data generating parameters of the simulation study using the isometric log-ratio transformation. 
- `01f_simulation_development_sb.Qmd`: Justification and visualization of the data generating parameters of the simulation study using the stick-breaking transformation.
- `02_fit_model_verbal_quantifiers.Qmd`: Empirical application of the Interval Truth Model to verbal quantifiers data, including model fitting, posterior analysis, and visualization of consensus intervals. The corresponding `.html` file shows a rendered version of the document.
	- Set `refit: true` in the YAML header if you want to refit all models upon rendering of the quarto markdown.
- `03_illustrations_article.Qmd`: Script that renders illutrating plots from the article.
- `04_prior_checks.Qmd`: Prior predictive checks and visualization of model priors to ensure reasonable model specification. The corresponding `.html` file shows a rendered version of the document.
- `05_zero_handling_sensitivity_checks.Qmd`: Sensitivity analysis for handling zero values in interval responses and comparison of different approaches. The corresponding `.html` file shows a rendered version of the document.


#### `src/simulation_v1_v2/`
Contains code and results for both versions of the simulation study:
- **v1**: First version of the simulation study with initial parameter settings
- **v2**: Updated simulation study with refined parameters and additional analyses
See the manuscript for more details on the differences between versions. 

#### `src/models`
Contains the Stan models.

- `itm_simulation_v2_beta.stan`: Stan model for the main simulation study (version 2).
- `itm_beta.stan`: Generic Stan model using the beta prior on the latent consensus intervals
- `itm_dirichlet.stan`: Generic Stan model using the Dirichlet prior on the latent consensus intervals
- `itm_quantifier_beta.stan`: Custom Stan model for the empirical dataset of verbal quantifiers


## Supplementary Materials Guide

This section provides a mapping between supplementary materials mentioned in the manuscript and their locations in this repository. Materials are organized by manuscript section for easy navigation.

### Introduction & Theory (Sections 1-2)

**Parameter Illustrations & Theory:**
- Parameter illustration figures → `plots/icm_parameter_illustration.pdf`
- Logit transformation illustration → `plots/logit_illustration.pdf`

**Prior Distributions & Sensitivity Analysis:**
- Prior predictive checks → `src/04_prior_checks.Qmd` and `src/04_prior_checks.html`
- Sensitivity analysis for zero-handling → `src/05_zero_handling_sensitivity_checks.Qmd` and `src/05_zero_handling_sensitivity_checks.html`
- Alternative Dirichlet prior implementation → `src/models/itm_dirichlet.stan`

### Simulation Study (Section 3)

**Main Simulation Study:**
- Simulation code → `src/01a_simulation_study_final_rev1.Rmd`
- Simulation results → `sim_results/sim_results_itm_2024-09-07_07-27-33/`
- Simulation visualizations → `src/01b_sim_visualizations_rev1.Qmd` and `src/01b_sim_visualizations_rev1.html`

**Link Function Study:**
- Link function simulation → `src/01c_simulation_study_link_functions_final.Rmd`
- Link function results → `sim_results/sim_res_link_function_2024-08-09_07-38-27/`
- Link function visualizations → `src/01d_sim_link_functions_visualizations.Qmd` and `src/01d_sim_link_functions_visualizations.html`

**Simulation Plots & Analysis:**
- Main simulation plots → `plots/sim_main/` (bias, MSE, scatterplots)
- Link function plots → `plots/sim_link_function/`
- Joint bias plots → `plots/sim_main/sim_absbias_widloc_comb.pdf`
- Raw estimation results → `plots/sim_main/raw_ests_boxplot.pdf`

**Parameter Development & Justification:**
- ILR transformation development → `src/01e_simulation_development_ilr.Qmd` and `src/01e_simulation_development_ilr.html`
- Stick-breaking transformation development → `src/01f_simulation_development_sb.Qmd` and `src/01f_simulation_development_sb.html`

### Empirical Example (Section 4)

**Main Empirical Analysis:**
- Full analysis code → `src/02_fit_model_verbal_quantifiers.Qmd`
- Rendered analysis report → `src/02_fit_model_verbal_quantifiers.html`
- Model fits → `fits/itm_quantifier_beta_custom_fit.RDS` and related files

**Empirical Plots & Results:**
- Consensus intervals plots → `plots/empirical_example/consensus_intervals_custom_model.pdf`
- Prior vs. posterior comparison → `plots/empirical_example/prior_vs_posterior_quantifiers_selection.pdf`
- All verbal quantifier densities → `plots/empirical_example/verbal_quantifier_densities_all.pdf`
- Example quantifier plots → `plots/empirical_example/verbal_quantifier_densities_example.pdf`
- Proficiency analysis → `plots/empirical_example/proficiency_vs_intervals.pdf`

**Sensitivity Analyses:**
- Zero-handling comparison → `plots/comparison_zero_handling_empirical.pdf`
- Consensus intervals sensitivity → `plots/consensus_intervals_sensitivity.pdf`

### Technical Implementation

**Stan Models:**
- Main simulation model → `src/models/itm_simulation_v2_beta.stan`
- Generic models → `src/models/itm_beta.stan`, `src/models/itm_dirichlet.stan`
- Empirical analysis model → `src/models/itm_quantifier_beta.stan`

**Session Information & Package Versions:**
- Session info provided in rendered HTML reports (`.html` files)
- Computational environment details included in analysis outputs



## Reproducing our results

If you want to reproduce our results, you first need to download the fit files 
and simulation results by running the script `src/00_download_big_files_osf.R`.
These files are too large to be synced via GitHub and therefore are stored in the 
OSF storage.


If you have installed Docker and Make, you can use the following files to reproduce the main results of our simulation study within a Docker container:

Run `make docker` from the root directory of this git repository. This will install all necessary
dependencies using the `Dockerfile` and `Makefile`. RStudio Server can then be opened from a browser
(<http://localhost:8787>), and the analysis documents can then be rerun.






