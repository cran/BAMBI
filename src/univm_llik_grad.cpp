#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif


#include "thread_num.h"
#include "bessel.h"
#include "univmgen.h"

typedef std::vector<arma::mat> stdvec_mat;


// [[Rcpp::export]]
double const_univm(double k) {
  return(2 * M_PI * BESSI0(k));
}

// [[Rcpp::export]]
arma::vec log_const_univm_all(arma::mat par_mat)
{
  int K = par_mat.n_cols;
  arma::vec result(K);
  for(int j = 0; j < K; j++)
    result[j] = log(const_univm(par_mat(0, j)));
  return(result);
}



// [[Rcpp::export]]
double ldunivmnum(double x, arma::vec par)
{
  return (par[0]*cos(x-par[1]));
}


// [[Rcpp::export]]
arma::vec dunivm_manyx_onepar(arma::vec x, double k, double mu)
{
  int n = x.n_rows;

  double l_const = log(const_univm(k));
  arma::vec par(2);
  par[0] = k; par[1] = mu;

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldunivmnum(x[i], par);
  }

  return arma::exp(ld_num - l_const);
}

// [[Rcpp::export]]
arma::vec dunivm_manyx_manypar(arma::vec x, arma::vec k, arma::vec mu)
{
  int n = k.size();

  arma::mat all_par(2, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k[i];
    all_par(1,i) = mu[i];
  }

  arma::vec l_const_all = log_const_univm_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldunivmnum(x[i], all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec dunivm_onex_manypar(double x, arma::vec k, arma::vec mu)
{
  int n = k.size();

  arma::mat all_par(2, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k[i];
    all_par(1,i) = mu[i];
  }

  arma::vec l_const_all = log_const_univm_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldunivmnum(x, all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}




//  [[Rcpp::export]]
arma::vec runivm_onepar(int n, double k, double mu)
{
  if(n == 1){
    arma::vec single_sample(1);
    single_sample[0] = runivm_single_onepar(k, mu);
    return(single_sample);
  } else {
    arma::vec all_sample(n);
    for(int i = 0; i < n; i++)
      all_sample[i] = runivm_single_onepar(k, mu);
    return(all_sample);
  }
}


//  [[Rcpp::export]]
arma::mat runivm_manypar(arma::vec k, arma::vec mu)
{
  int n = k.size();
  arma::vec all_sample(n);
  for(int i = 0; i < n; i++)
    all_sample[i] = runivm_single_onepar(k[i], mu[i]);
  return(all_sample);
}


// [[Rcpp::export]]
double univmmix(double x,  arma::mat par, arma::vec pi, arma::vec log_c_von)
{
  double res = 0;
  int K = par.n_cols;

  for(int j = 0; j < K; j++)
    res += exp(ldunivmnum(x, par.col(j)) - log_c_von[j]) * pi[j];
  return(res);
} //likelihood contribution (density) of a single point in a mixture model


// [[Rcpp::export]]
arma::vec univmmix_manyx(arma::vec x, arma::mat par, arma::vec pi, arma::vec log_c)
{
  int n = x.n_rows;
  arma::vec result(n);
  for(int i = 0; i < n; i++)
    result[i] = univmmix(x[i], par, pi, log_c);
  return(result);
}


// [[Rcpp::export]]
arma::mat mem_p_univm(arma::vec data, arma::mat par, arma::vec pi, arma::vec log_c_von, int ncores = 1)
{
  int n = data.n_rows, K = par.n_cols, j;
  double row_total;
  arma::mat den(n, K);
#pragma omp parallel for private(j, row_total) num_threads(ncores)
  for(int i = 0; i < n; i++){
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = pi[j]*exp(ldunivmnum(data[i], par.col(j)) - log_c_von[j]);;
      row_total += den(i, j);
    }
    row_total = maxi(row_total, 1e-50);
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);

}


// [[Rcpp::export]]
double llik_univm_full(arma::vec data, arma::mat par, arma::vec pi,
                       arma::vec log_c, int ncores = 1)
{
  int n = data.n_rows, K = pi.size(), j;
  long double temp, log_sum = 0.0;
  arma::vec log_pi = log(pi);

  if(K > 1) {
#pragma omp parallel for reduction(+:log_sum)  private(j, temp) num_threads(ncores)
    for(int i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(ldunivmnum(data[i], par.col(j)) - log_c[j] + log_pi[j] );
      log_sum += log(maxi(temp, 1e-100));
    }
  } else {
#pragma omp parallel for reduction(+:log_sum)  num_threads(ncores)
    for(int i = 0; i < n; i++)
      log_sum += ldunivmnum(data[i], par);
    log_sum -= n*log_c[0];

  }
  return(log_sum);
}


arma::vec grad_log_univm_one_comp_i_unadj(double x, arma::vec par, double bes_ratio_1_0_kappa)
{
  double mu = par[1];
  arma::vec result(2);
  result[0] = cos(x - mu) - bes_ratio_1_0_kappa;
  result[1] = sin(x - mu);
  return (result);
} // without the multiplier kappa in del_mu



// [[Rcpp::export]]
arma::mat grad_univm_all_comp(arma::vec data, arma::mat par_mat, arma::vec pi, int ncores = 1)
{
  int n = data.size(), K = pi.size(), j;
  double denom, temp;
  arma::mat grad_temp(2, K), grad_sum = arma::zeros(2, K);

  arma::vec bes_1_kappa(K), bes_0_kappa(K), bes_ratio_1_0_kappa(K), log_pi = log(pi);

  for(j = 0; j < K; j++) {
    bes_0_kappa[j] = BESSI0(par_mat(0, j));
    bes_1_kappa[j] = BESSI1(par_mat(0, j));
    bes_ratio_1_0_kappa[j] = bes_1_kappa[j]/bes_0_kappa[j];
  }

  arma::vec log_univm_const = log(2 * M_PI * bes_0_kappa);

  stdvec_mat grad_sum_part(ncores);
  for(int g = 0; g < ncores; g++)
    grad_sum_part[g] = arma::zeros(2, K);


#pragma omp parallel for   private(j, denom, temp) num_threads(ncores)
  for(int i = 0; i < n; i++) {
    arma::mat grad_temp(2, K);
    int g = get_thread_num_final();
    denom = 0;
    for(j = 0; j < K; j++) {
      temp = exp(ldunivmnum(data[i], par_mat.col(j)) - log_univm_const[j] + log_pi[j] );
      denom += temp;
      grad_temp.col(j) = temp*grad_log_univm_one_comp_i_unadj(data[i], par_mat.col(j), bes_ratio_1_0_kappa[j]);
    }
    grad_sum_part[g] += grad_temp/denom;
  }  // the pi_h term in numerator hasn't been adjusted yet!

  for(int g = 0; g < ncores; g++)
    grad_sum += grad_sum_part[g];

  arma::mat adj_factor = arma::ones(2, K);
  adj_factor.row(1) = par_mat.row(0);

  return (grad_sum % adj_factor);
}



