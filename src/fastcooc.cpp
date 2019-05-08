// 
// Build a fast coocurrence network
// 

#include <RcppArmadillo.h>
using namespace arma; 

//[[Rcpp::export]]
Rcpp::List cooc_mat(const arma::umat comm, 
                    const arma::cube pvalarr) { 
  
  // const uword nrep = 99; 
  uword nsp = comm.n_cols; 
  
  // arma::umat coocmat(nsp, nsp, fill::zeros); 
  arma::mat pmat(nsp, nsp, fill::zeros);  
  
  for (uword spi=0; spi<nsp; spi++) { 
    
    // Set diag to non-significant value
    pmat(spi, spi) = 0.5; 
    
    for (uword spj=(spi+1); spj<nsp; spj++) { 
      
      uword total_i = accu(comm.col(spi)); 
      uword total_j = accu(comm.col(spj)); 
      
      uword coocobs = accu( comm.col(spi) % comm.col(spj) ); 
      
      // coocmat(spi, spj) = coocobs; 
      // Rcpp::Rcout << total_i << " " << total_j << " " << coocobs << " \n"; 
      
      pmat(spi, spj) = pvalarr(total_i, 
                               total_j, 
                               coocobs); // here coocobs is between 0 and 12
                                         // so we can use it directly as an 
                                         // index
      pmat(spj, spi) = pmat(spi, spj); 
      
      // arma::uvec cooctries(nrep+1); 
      // cooctries(1) = obs; 
      // for ( uword n=1; n<(nrep+1); n++ ) { 
      //   cooctries(n) = accu( shuffle(comm.col(spi)) % 
      //                        shuffle(comm.col(spj)) ); 
      // }
      // This return where in the vector we should put ix 
      // so that the vector is sorted. In other words 
      // it is its rank. 
      // uvec ix = sort_index(cooctries); 
      // pmat(spi, spj) = ix(1); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("pvalues") = pmat); 
}

//[[Rcpp::export]]
arma::vec coocvec(const arma::uvec veci, 
                  const arma::uvec vecj, 
                  const arma::uword ntries
                 ) { 
  
  const uword total_i = accu(veci); 
  const uword total_j = accu(vecj); 
  const uword coocobs = accu(veci % vecj); 
  
  double below = 0; 
  for (uword i=0; i<ntries; i++) { 
    uword cooctry = accu( shuffle(veci) % shuffle(vecj) ); 
    // Rcpp::Rcout << "try: " << cooctry << " / obs: " << coocobs << "\n"; 
    below += cooctry < coocobs ? 1 : 0; 
    below += cooctry == coocobs ? 0.5 : 0; 
  }
  
  double pval = double(below) / (1 + double(ntries)); 
  return arma::vec { double(total_i), 
                     double(total_j), 
                     double(coocobs), 
                     pval }; 
}
