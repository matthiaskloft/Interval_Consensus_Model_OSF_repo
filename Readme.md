# Interval Truth Model

This repository contains code and supplementary files for the manuscript:

> Kloft, M.,  Siepe, B.S., Heck, D.W. (2024). The Interval Truth Model: A Consensus Model for Continuous Bounded Interval Responses. <TODO ADD LINK>.

A BiBTeX entry for LaTeX users is:

```BibTeX
@article{Kloft2024,
  title={The Interval Truth Model: A Consensus Model for Continuous Bounded Interval Responses},
  author={Kloft, Matthias and Siepe, Bj{\"o}rn S. and Heck, Daniel W.},
  url={TODO},
  year={2024},
  note={PsyArXiv preprint}
}
``` 




## Folder Structure
### `data/`

This folder contains the data used in the empirical example.


### `plots/`

This folder contains all plots used in the manuscript and in the supplementary material.

### `sim_results/`

This folder contains the results of the simulation studies. 

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
- `02_fit_model_verbal_quantifiers.Qmd`: TODO. The corresponding `.html` file shows a rendered version of the document.
	- Set `refit: true` in the YAML header if you want to refit all models upon rendering of the quarto markdown.
- `03_illustrations_article.Qmd`: Script that renders illutrating plots from the article.
- `04_prior_checks.Qmd`: Visualization of XXX. The corresponding `.html` file shows a rendered version of the document.


#### `src/simulation_v1`
Contains code and results for the first version of the simulation study. See the manuscript for more details. 

#### `src/models`
Contains the Stan models.

- `itm_simulation_v2_beta.stan`: Stan model for the main simulation study (version 2).
- ``: Generic Stan model using the beta prior on the latent consensus intervals
- ``: Generic Stan model using the Dirichlet prior on the latent consensus intervals
- ``: Custom Stan model for the empirical dataset of verbal quantifiers


## Reproducing our results

If you want to reproduce our results, you first need to download the fit files 
and simulation results by running the script `src/00_download_big_files_osf.R`.
These files are too large to be synced via GitHub and therefore are stored in the 
OSF storage.


If you have installed Docker and Make, you can use the following files to reproduce the main results of our simulation study within a Docker container:

Run `make docker` from the root directory of this git repository. This will install all necessary
dependencies using the `Dockerfile` and `Makefile`. RStudio Server can then be opened from a browser
(<http://localhost:8787>), and `05_sim_visualizations_v2.Qmd` can then be rerun.


## Note
In some places of the project, you might encounter the acronym ITM/itm, which stands for Interval Truth Model. 
This is the former name of the model, which we changed to Interval Consensus Model (ICM) in the final version of the manuscript.



