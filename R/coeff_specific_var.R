#' @title coeff_specific_var
#'
#' @description A function to retrieve a vector of coefficient estimates for a specific variable of interest after running PRANA.
#'
#' @param  coefftab A table that includes adjusted p-values for a specific variable.
#' @param varname Specify the name of the variable of interest.
#'
#' @return A vector of coefficient estimates for a single variable from the model.
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
#' # Create an object to keep the table with coefficient estimates using coeff() function.
#' coefftab <- coeff(PRANAres)
#'
#' # Lastly, we use coeff_specific_var function to retrieve
#' # adjusted p-values for a single variable of interest.
#' coeff_specific_var(coefftab = coefftab, varname = "currentsmoking")
#' @export

coeff_specific_var <- function(coefftab, varname) {
        coeffvec <- coefftab[varname]
        return(coeffvec = coeffvec)
}
