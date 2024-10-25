packages <- c(
  "here",
  "osfr"
)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(packages, update = F, character.only = T)


# remove existing fits folders
if (dir.exists(here("fits"))) {
  unlink(here("fits"), recursive = T)
}
if (dir.exists(here("sim_results"))) {
  unlink(here("sim_res"), recursive = T)
}


# Specify OSF files
project <- osf_retrieve_node("https://osf.io/r32by")
files <- osf_ls_files(project)
# download
osf_download(files[], )

