## set R version (https://hub.docker.com/r/rocker/verse/tags)
FROM rocker/verse:4.4.1

## set up directories
RUN mkdir /home/rstudio/sim_results /home/rstudio/src 
COPY interval_truth_model.Rproj /home/rstudio/

## install R packages from CRAN the last day of the specified R version
## ncpus set to -1 (all available cores)
RUN install2.r --error --skipinstalled --ncpus -1 \
    tidyverse SimDesign rstan posterior bayesplot here pacman psych ggh4x ggokabeito ggExtra showtext ggdist pander sysfonts knitr

