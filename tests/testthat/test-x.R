# Load and process the data to run PRANA
data(combinedCOPDdat_RGO)
rnaseqdat = combinedCOPDdat_RGO[ , 8:ncol(combinedCOPDdat_RGO)]
rnaseqdat = as.data.frame(apply(rnaseqdat, 2, as.numeric))

phenodat = combinedCOPDdat_RGO[order(combinedCOPDdat_RGO$currentsmoking), 2:7]

index_grpA = which(combinedCOPDdat_RGO$currentsmoking == 0)
index_grpB = which(combinedCOPDdat_RGO$currentsmoking == 1)

PRANAres <- PRANA(RNASeqdat = rnaseqdat, clindat = phenodat,
                  groupA = index_grpA, groupB = index_grpB)

# Create a function to check the sign of p-values
# and adjusted p-values:
sign_check <- function(pval) {
        if (sum(pval < 0) > 0) {
                return("ERROR: Negative p-values")
        } else {
                return("No negative p-values")
        }
}


# A function to check if results have any negative p-values:
test_that("Check any negative p-values after PRANA", {
        expect_equal(sign_check(PRANAres$p_values), "No negative p-values")
        expect_equal(sign_check(PRANAres$adjp_values), "No negative p-values")
})

# Use testthat to check if results return any missing or NA:
# Note: coefficient estimates and p-values
#       (and adjusted p-values) should not have any missing or NA.
test_that("Check the missingness after PRANA", {
        expect_false(any(is.na(PRANAres$beta_hat)))
        expect_false(any(is.na(PRANAres$p_values)))
        expect_false(any(is.na(PRANAres$adjp_values)))
})









