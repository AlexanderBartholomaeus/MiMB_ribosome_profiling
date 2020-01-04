### script to install R packages

# install CRAN packages
packs <- c(
  'foreach',
  'doParallel'
)
install.packages(packs)

# install Bioconductor packages
