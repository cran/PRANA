#' @title coeff
#'
#' @description A function to retrieve a table with coefficient estimates after running PRANA.
#'  The table includes all variables that were included in the pseudo-value regression model.
#'
#' @param PRANAres An object called after running PRANA.
#'
#' @return A table that includes coefficient estimates for all variables included in the fitted model.
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
#' # Then, create an object called `res` to call results later.
#' PRANAres <- PRANA(RNASeqdat = rnaseqdat, clindat = phenodat,
#'  groupA = index_grpA, groupB = index_grpB)
#'
#' # Now, we want to keep the table with coefficient estimates only.
#' coeff(PRANAres)
#' @export

coeff <- function(PRANAres) {
        coefftab <- PRANAres$beta_hat
        return(coefftab = coefftab)
}
