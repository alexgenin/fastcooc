# 
# This function creates an array containing the P-values of observing at least 
# a certain number of co-occurrences given the total number of 
# occurrence of two species in a pair. 
# 
precompute_pvalues <- function(n, ntries = 1999, memoise = TRUE) { 
  
  pval_arr <- array(NA_real_, dim = c(n+1, n+1, n+1), 
                    dimnames = list(total_i = seq(0, n), 
                                    total_j = seq(0, n), 
                                    obscooc = seq(0, n)))
  for ( max_total in seq(0, n) ) { 
    for ( min_total in seq(0, max_total) ) { 
      for ( obscooc in seq(0, min_total) ) { 
        # If one of the two species is always present or always absent, we 
        # cannot say anything about how one distributes itself in space 
        # compared to the other. For those cases, we set the P-value to 0.5
        # (non-significant) in the 'else' statement. 
        if ( min_total != 0 && 
             max_total != 0 && 
             min_total != n && 
             max_total != n ) { 
          
          # We build a vector of presence/absence that respects the total 
          # abundance of the less abundant species and the more abundant 
          # species, as well as their number of co-occurrences. We then 
          # use that as a vector of presence/absence from which we compute 
          # 
          if ( (n - max_total) - min_total + obscooc >= 0 ) { 
            a <- c(rep(TRUE, max_total), rep(FALSE, n - max_total))
            b <- c(rep(TRUE,  obscooc), 
                  rep(FALSE, max_total - obscooc), 
                  rep(TRUE,  min_total - obscooc), 
                  rep(FALSE, (n - max_total) - min_total + obscooc))
            pval <- coocvec(a, b, ntries)[4]
            pval_arr[min_total+1, max_total+1, obscooc+1] <- pval
            pval_arr[max_total+1, min_total+1, obscooc+1] <- pval
          }
          
          # coocvec(a, b)
          
        } else { 
          pval_arr[min_total+1, max_total+1, obscooc+1] <- 0.5
          pval_arr[max_total+1, min_total+1, obscooc+1] <- 0.5
        }
      }
    }
  }
  return(pval_arr)
}
