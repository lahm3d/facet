# install_packages: This is code to installs necessary 
# packages before executing
# post processing scripts in R.

# define CRAN repository to download packages
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

# list of packages used for post-processing
list_of_packages <- c("dplyr", "foreign", "pastecs")

# get a list of all packages and check if installed or not installed
packages_to_install <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

# install relevant packages
if(length(packages_to_install)) install.packages(packages_to_install)