install.packages('PRROC', repos="https://cran.r-project.org/")
install.packages('data.table', repos="https://cran.r-project.org/")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cran.r-project.org/")
BiocManager::install(c("slingshot", "MAST"))
