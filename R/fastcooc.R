# 
# Wrapper around cooc_mat to put results in the right format
# 
# Comm_matrix is the P/A community matrix (sites as rows)
# Pval_array is the precomputed array of P-values
# long_form controls whether to return a (dense) matrix or a data.frame
# 
fastcooc <- function(comm_matrix, 
                     pval_array = NULL, 
                     long_form = TRUE, 
                     trim_threshold = 0.05, 
                     ...) { 
  
  if ( !is.logical(comm_matrix) ) { 
    stop("comm_matrix must be logical (TRUE/FALSE values")
  }
  
  if ( is.null(pval_array) ) { 
    message('P-value lookup array has not been pre-computed, computing it now')
    pval_array <- precompute_pvalues(nrow(comm_matrix), ...)
  }
  
  if ( nrow(comm_matrix) != (nrow(pval_array) - 1) || 
       nrow(comm_matrix) != (ncol(pval_array) - 1 )) { 
    stop('Array dimension mismatch. Check that you have ', nrow(pval_array), 
         'sites (rows) in your community matrix')
  }
  
  if ( is.null(colnames(comm_matrix)) ) { 
    spnames <- seq.int(ncol(comm_matrix))
  } else { 
    spnames <- colnames(comm_matrix)
  }
  
  cooc_results <- cooc_mat(comm_matrix, pval_array)$pvalues
  rownames(cooc_results) <- colnames(cooc_results) <- spnames
  
  if ( long_form ) { 
    cooc_results <- data.frame(expand.grid(sp1 = spnames, 
                                           sp2 = spnames), 
                               pval = as.vector(cooc_results))
    if ( trim_threshold < 0.5 ) { 
      cooc_results <- subset(cooc_results, 
                             pmin(pval, 1 - pval) < (trim_threshold/2))
    }
  }
  
  return(cooc_results)
}
# 
# Nsp <- 9000
# bigmat <- matrix(rnorm(12*Nsp), nrow = 12) > 0
# bigmat <- bigmat[ ,apply(bigmat, 2, sum) > 0] 
