# prerequisite packages
requiredPackages.cran <- c(
  "reshape2",
  "tidyr",
  "ggplot2")

requiredPackages.bioconductor <- c(
  "ComplexHeatmap")

# install packages if not installed on your system
# CRAN
install.packages(requiredPackages.cran)
# bioconductor 
source("https://bioconductor.org/biocLite.R")
biocLite(requiredPackages.bioconductor)
