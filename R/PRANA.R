#' @title PRANA
#'
#' @description A pseudo-value regression approach for differential network analysis
#'     that adjusts for additional covariates (PRANA)
#'
#' @param RNASeqdat An RNA-Seq data with subjects in rows and genes in columns.
#' @param clindat A data with clinical variables to be included in the regression
#'               (e.g., binary group variable indicating
#'               current smoking status, continuous age, ...)
#' @param groupA Indices of the subjects in the first category (e.g., non-current smoker)
#'               of binary group variable.
#' @param groupB Indices of the subjects in the second category (e.g., current smoker)
#'               of binary group variable.
#'
#' @return A list containing three data frame objects that summarize the results of PRANA. This includes
#'         beta coefficients, p-values, and adjusted p-values via the empirical Bayes approach
#'         for each predictor variables that are included in the regression model.
#'
#' @references Ahn S, Grimes T, Datta S. A pseudo-value regression approach for differential network analysis of co-expression data. BMC Bioinformatics, 2023;24(1):8
#'
#' @examples
#' data(combinedCOPDdat_RGO) # A complete data containing expression and clinical data.
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
#' PRANAres <- PRANA(RNASeqdat = rnaseqdat, clindat = phenodat,
#' groupA = index_grpA, groupB = index_grpB)
#' @export
#' @import dnapath
#' @import dplyr
#' @import parallel
#' @import robustbase
#' @import minet

PRANA <- function(RNASeqdat, clindat, groupA, groupB) {
        est_method = run_aracne # ARACNE are used to estimate the network (association matrix)

        #############################################################################################
        #  STEP 1. Estimate an association matrix via ARACNE from the RNA-seq expression data.
        #############################################################################################
        rnaseqdatA = RNASeqdat[groupA, ] # RNA-Seq data for non-current smoker (namely Group A)
        rnaseqdatB = RNASeqdat[groupB, ] # RNA-Seq data for current smoker (namely Group B)
        nw_est_grpA = est_method(rnaseqdatA, verbose = F) # Association matrix for Group A
        nw_est_grpB = est_method(rnaseqdatB, verbose = F) # Association matrix for Group B

        n_A <- length(groupA) # Sample size for Group A
        n_B <- length(groupB) # Sample size for Group B


        #############################################################################################
        #  STEP 2. Calculate \hat{\theta}_{k} (total connectivity) by taking the column sum of the
        #          association matrix to obtain the total connectivity of for each gene.
        #############################################################################################
        thetahat_grpA = thetahats(nw_est_grpA)
        thetahat_grpB = thetahats(nw_est_grpB)


        #############################################################################################
        #  STEP 3. Re-estimate association matrix using the expression data without i-th subject.
        #          Then, calculate \hat{\theta}_{k(i)} from the reestimated association matrix.
        #############################################################################################
        # Re-estimation part
        nw_est_drop_grpA <- mclapply(groupA, function(j) est_method(rnaseqdatA[-j, ], verbose = F))
        nw_est_drop_grpB <- mclapply(groupB, function(j) est_method(rnaseqdatB[-j, ], verbose = F))
        # thetahat_{-i} for each gene
        thetahat_drop_grpA <- sapply(nw_est_drop_grpA, thetahats)
        thetahat_drop_grpB <- sapply(nw_est_drop_grpB, thetahats)



        #############################################################################################
        #  STEP 4. Calculate jackknife pseudovalues (\tilde{\theta}_{ik} using \hat{\theta}_{k}
        #          and \hat{\theta}_{k(i)}.
        #############################################################################################
        thetatildefun <- function(thetahatinput, thetahatdropinput, sizegroup) {
                thetatildeout = matrix(NA, ncol=length(thetahatinput), nrow=sizegroup)
                thetatildeout = sapply(1:nrow(thetahatdropinput), function(k) {
                        sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
                })
                return(thetatildeout)
        }

        thetatilde_grpA = thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
        thetatilde_grpB = thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)
        thetatilde = rbind(thetatilde_grpA, thetatilde_grpB)
        colnames(thetatilde) = colnames(RNASeqdat) # Map the column names (gene names)


        #############################################################################################
        #  STEP 5. Fit a robust regression model
        ##############################################################################################
        pseudo.beta_list <- lapply(1:ncol(thetatilde), function(i) {
                m <- thetatilde[, i]
                df <- data.frame(clindat,
                                 m = m)
                fit <- ltsReg(m ~currentsmoking + packyrs + age + gender + race + FEV1perc, data = df, mcd=FALSE) # Include a set of covariates to be regressed.
                return(fit)
        })

        ### Obtain p-values for each genes:
        beta_hat = vector(mode = "list", ncol(thetatilde))
        p_values = vector(mode = "list", ncol(thetatilde))
        k = NULL
        for(k in 1:ncol(thetatilde)) {
                # Estimates for the beta coefficients:
                beta_hat[[k]] <- summary(pseudo.beta_list[[k]])$coef[-1, "Estimate"]
                p_values[[k]] <- summary(pseudo.beta_list[[k]])$coef[-1, "Pr(>|t|)"]

        }

        beta_hat = as.data.frame(bind_rows(beta_hat))  # Convert list into data.frame
        rownames(beta_hat) <- colnames(RNASeqdat) # Map the gene names to the data.frame for betahats

        p_values = as.data.frame(bind_rows(p_values))  # Convert list into data.frame
        rownames(p_values) <- colnames(RNASeqdat) # Map the gene names to the data.frame for p-values

        # Compute the adjusted p-values via empirical Bayes approach proposed by the reference below:
        # NOTE: EBS() is code from Datta S and Datta S (2005). Empirical Bayes screening of many p-values with applications to microarray studies. Bioinformatics, 21(9), 1987-1994
        adjp_values = as.data.frame(apply(p_values, 2, function(x) EBS(pvo = x, alpha = 0.05, B = 500, h = 1)))
        rownames(adjp_values) <- colnames(RNASeqdat) # Map the gene names to the data.frame for adj p-values

        return(list(beta_hat = beta_hat, p_values = p_values, adjp_values = adjp_values))
}
