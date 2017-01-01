#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec rowVars(arma::mat mat_in)
{
  int nrow = mat_in.n_rows;
  arma::vec res(nrow);
  for(int i = 0; i < nrow; i++) {
    res[i] = arma::var(mat_in.row(i));
  }
  return res;
}


// [[Rcpp::export]]
arma::cube par_mat_permute(arma::cube par_mat, arma::umat perm_lab)
{
  int n_iter = par_mat.n_slices, n_row = par_mat.n_rows, n_col = par_mat.n_cols;
  arma::cube result(n_row, n_col, n_iter);
  for(int iter = 0; iter < n_iter; iter++) {
    for(int row = 0; row < n_row; row++) {
      for(int col = 0; col < n_col; col++) {
        result(row, col, iter) = par_mat(row, perm_lab(iter, col) - 1, iter);
      }
    }
  }

  return result;
}


// [[Rcpp::export]]
Rcpp::NumericVector cID(arma::mat probs, int ncomp, arma::vec Uv) {
  double U;
  double* p = new double[ncomp];
  Rcpp::NumericVector clID(probs.n_rows);

  for (int i = 0; i < (int)probs.n_rows; i++) {
    U = Uv[i];
    //Rcout << U;
    p[0] = probs(i,0);
    if (U < p[0]) {
      clID[i] = 1;
      continue;
    }

    for (int j = 1; j < ncomp; j++) {
      p[j] = probs(i,j) + p[j-1];
      if (U < p[j]) {
        clID[i] = j+1;
        break;
      }
    }

  }

  delete[] p;
  return(clID);
}
