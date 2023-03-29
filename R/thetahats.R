#' @title thetahats
#'
#' @description A function to compute the total connectivity of each gene from the association matrix.
#'
#' @param asso.matinput An association matrix that is estimated from the user-provided expression data is used
#'                      as an input to compute the total connectivity of each gene.
#'
#' @return A vector containing total connectivity of each gene (i.e. continuous version of centrality measure of a network)
#'
#' @export

thetahats = function(asso.matinput) {
        results = vector()
        for(j in 1:ncol(asso.matinput)) {
                results[j] = sum(asso.matinput[j, -j])
        }
        return(results)
}
