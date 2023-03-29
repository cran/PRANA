#' @title sigDCGtab
#'
#' @description A function to retrieve the data frame that are significantly differentially connected (DC).
#' between two biological/clinical states (aka the main binary indicator) with the presence of additional covariate information.
#'
#' @param adjptab A table with adjusted p-values and names for the variable that the user specifies in the groupvar.
#' @param groupvar Specify the name of binary indicator variable.
#' @param alpha A level of significance (e.g. 0.05).
#'
#' @return Adjusted p-values and names of significantly DC genes (e.g. gene IDs) from PRANA.
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
#'  groupA = index_grpA, groupB = index_grpB)
#'
#' # Next, we want to call the table with adjusted p-values.
#' adjptab <- adjpval(PRANAres)
#'
#' # Please specify the name of variable in sigDCGtab(groupvar = ).
#' sigDCGtab <- sigDCGtab(adjptab = adjptab, groupvar = "currentsmoking", alpha = 0.05)
#' sigDCGtab
#' @export

sigDCGtab <- function(adjptab, groupvar, alpha) {
        adjpvalvec = adjptab[groupvar]
        sigDCGtab = filter(adjpvalvec, adjpvalvec < alpha)
        return(sigDCGtab = sigDCGtab)
}

