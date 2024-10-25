packages <- c("here", "osfr")

if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(packages, update = F, character.only = T)


# Specify OSF file
project <- osf_retrieve_node("https://osf.io/r32by")
file <- osf_ls_files(project, verbose = TRUE)

# download
osf_download(file, conflicts = "overwrite", verbose = TRUE)
# unzip packaged files
unzip(here("fits_and_sim_results.zip"), overwrite = TRUE)
unlink(here("fits_and_sim_results.zip"))
