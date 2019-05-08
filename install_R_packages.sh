#!/usr/bin/env bash

module load gcc
module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

# Install data.table in R version 3.5.1
if (!require(data.table)) install.packages("data.table")

# Install data.table in R version 3.5.1
if (!require(foreach)) install.packages("foreach")

# Install doMC in R version 3.5.1
if (!require(doMC))  install.packages("doMC")

# Install ggplot2 in R version 3.5.1
if (!require(ggplot2))  install.packages("ggplot2")

# Check BiocManager is installed in R version 3.5.1
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install GWASTools in R Version 3.5.1
if (!require(GWASTools)) BiocManager::install("GWASTools", version = "3.8")

# Install SNPRelate in R Version 3.5.1
if (!require(SNPRelate)) BiocManager::install("SNPRelate", version = "3.8")

# Install gdsfmt in R version 3.5.1
if (!require(gdsfmt)) BiocManager::install("gdsfmt", version = "3.8")

# Install GENESIS in R version 3.5.1
if (!require(GENESIS)) BiocManager::install("GENESIS", version = "3.8")

# Install apeglm in R version 3.5.1
if (!require(apeglm)) BiocManager::install("apeglm", version = "3.8")



library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
library(GWASTools)
library(SNPRelate)
library(gdsfmt)
library(GENESIS)
