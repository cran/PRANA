## ---- include = FALSE---------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# library(dplyr) 
# library(parallel) # To use mclapply() when reestimating the association matrix.
# library(robustbase) # To fit a robust regression
# library(minet) # To estimate network via ARACNE

# Please run the following lines to install minet package 
# from Bioconductor in your R console:
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("minet")

## ----setup--------------------------------------------------------------------
library(PRANA)

## -----------------------------------------------------------------------------
data(combinedCOPDdat_RGO) # A complete data containing expression and clinical data.

## -----------------------------------------------------------------------------
# A gene expression data part of the downloaded data.
rnaseqdat <- combinedCOPDdat_RGO[ , 8:ncol(combinedCOPDdat_RGO)]
rnaseqdat <- as.data.frame(apply(rnaseqdat, 2, as.numeric))

# A clinical data with additional covariates sorted by current smoking groups:
# The first column is ID, so do not include.
phenodat <- combinedCOPDdat_RGO[order(combinedCOPDdat_RGO$currentsmoking), 2:7]

# Indices of non-current smoker (namely Group A)
index_grpA <- which(combinedCOPDdat_RGO$currentsmoking == 0)
# Indices of current smoker (namely Group B)
index_grpB <- which(combinedCOPDdat_RGO$currentsmoking == 1)

## -----------------------------------------------------------------------------
PRANAres <- PRANA(RNASeqdat = rnaseqdat, clindat = phenodat, groupA = index_grpA, groupB = index_grpB)

## -----------------------------------------------------------------------------
# This is useful when we want to create a table with adjusted p-values only.
adjpval(PRANAres)

## -----------------------------------------------------------------------------
# Create an object to keep the table with adjusted p-values using adjpval() function.
adjptab <- adjpval(PRANAres)

## -----------------------------------------------------------------------------
# NOTE: Please do NOT forget to provide a name of variable with the quotation marks!
adjpval_specific_var(adjptab = adjptab, varname = "currentsmoking")

## -----------------------------------------------------------------------------
# NOTE: Please do NOT forget to provide a name of variable with the quotation marks!
sigDCGtab <- sigDCGtab(adjptab = adjptab, groupvar = "currentsmoking", alpha = 0.05)
sigDCGtab

## -----------------------------------------------------------------------------
# NOTE: Please do NOT forget to provide a name of variable with the quotation marks!
sigDCGnames <- sigDCGnames(adjptab = adjptab, groupvar = "currentsmoking", alpha = 0.05)
sigDCGnames

## -----------------------------------------------------------------------------
#rename_genes(sigDCGnames, to = "symbol", species = "human")

