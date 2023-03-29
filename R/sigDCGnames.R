#' @title sigDCGnames
#'
#' @description A function to retrieve the name of genes that are significantly differentially connected (DC).
#' between two biological/clinical states (aka the main binary indicator) with the presence of additional covariate information.
#'
#' @param adjptab A table with adjusted p-values for all variables that were included in the pseudo-value regression model.
#' @param groupvar Specify the name of binary indicator variable.
#' @param alpha A level of significance (e.g. 0.05).
#'
#' @return Names of significantly DC genes (e.g. gene IDs) from PRANA.
#' If you need both adjusted p-values and names, please use sigDCGtab() instead.
#'
#' @examples
#' #' data(combinedCOPDdat_RGO) # A complete data containing expression and clinical data.
#'
#' # A gene expression data part of the downloaded data.
#' rnaseqdat = combinedCOPDdat_RGO[ , 8:ncol(combinedCOPDdat_RGO)]
#' rnaseqdat = as.data.frame(apply(rnaseqdat, 2, as.numeric))
#'
#' # A clinical data with additional covariates sorted by current smoking groups:
#' # The first column is ID, so do not include.
#' phenodat = combinedCOPDdat_RGO[order(combinedCOPDdat_RGO$currentsmoking), 2:7]
#'
#' # Indices of non-current smoker (namely Group A)
#' index_grpA = which(combinedCOPDdat_RGO$currentsmoking == 0)
#' # Indices of current smoker (namely Group B)
#' index_grpB = which(combinedCOPDdat_RGO$currentsmoking == 1)
#'
#' # Use PRANA() function to perform the pseudo-value regression analysis.
#' # Then, create an object called PRANA_Results to call results.
#' PRANAres <- PRANA(RNASeqdat = rnaseqdat, clindat = phenodat,
#' groupA = index_grpA, groupB = index_grpB)
#'
#' # Next, we want to call the table with adjusted p-values only.
#' adjptab <- adjpval(PRANAres)
#'
#' # Please specify the name of binary group indicator in sigDCGnames(groupvar = ).
#' sigDCGnames <- sigDCGnames(adjptab = adjptab, groupvar = "currentsmoking", alpha = 0.05)
#' sigDCGnames
#' @export

sigDCGnames <- function(adjptab, groupvar, alpha) {
        adjpvalvec = adjptab[groupvar]
        DCGnames = rownames(adjpvalvec)[adjpvalvec < alpha]
        return(DCGnames = DCGnames)
}
