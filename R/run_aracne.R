#' @title run_aracne
#'
#' @description A function to conduct co-expression analysis using ARACNE (Margolin et al. 2006). Uses the implementation from the minet package (Meyer et al. 2008).
#'     This function is from dnapath R package, which is archived in May 2024.
#'
#' @param x A n by p matrix of gene expression data (n samples and p genes)
#' @param estimator Argument is passed into minet::build.mim.
#' @param disc Argument is passed into minet::build.mim.
#' @param nbins Argument is passed into minet::build.mim.
#' @param eps Argument is passed into minet::aracne.
#' @param ... Additional arguments are ignored.
#'
#' @return A p by p matrix of association scores.
#'
#' @references Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A (2006). “ARACNE: An Algorithm for the Reconstruction of Gene Regulatory Networks in a Mammalian Cellular Context.” In BMC Bioinformatics, volume 7(1), S7. BioMed Central.
#' @references Meyer PE, Lafitte F, Bontempi G (2008). “minet: A R/Bioconductor Package for Inferring Large Transcriptional Networks using Mutual Information.” BMC Bioinformatics, 9(1), 461.
#'
#' @export
#' @import minet

run_aracne <- function (x, estimator = "spearman", disc = "none", nbins = NULL,
          eps = 0, ...)
{
        #if (!requireNamespace("minet", quietly = TRUE)) {
        #        warning("The `minet` package must be installed to use run_aracne(). Using ",
        #                "run_pcor() instead.")
        #        return(run_pcor(x, ...))
        #}
        p <- ncol(x)
        scores <- matrix(0, nrow = p, ncol = p)
        index <- which(apply(x, 2, function(val) !all(abs(val - mean(val)) <
                                                              1e-12)))
        if (length(index) <= 1)
                return(scores)
        mim <- minet::build.mim(x[, index], estimator = estimator,
                                disc = disc, nbins = nbins)
        scores[index, index] <- minet::aracne(mim, eps = eps)
        colnames(scores) <- colnames(x)
        rownames(scores) <- NULL
        return(scores)
}
