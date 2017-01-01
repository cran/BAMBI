#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include "bessel.h"
#include "univmgen.h"

typedef std::vector<arma::mat> stdvec_mat;
typedef std::vector<arma::vec> stdvec_vec;

// ---------------------------------------------------------------------
// VON MISES SINE
// ---------------------------------------------------------------------

// [[Rcpp::export]]
double const_vmsin(double k1, double k2, double lambda)
{

  double start = BESSI(0,k1) * BESSI(0,k2), rat = pow(lambda, 2) / (4*k1*k2);
  int p = 1;
  double temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
  start += temp;
  while((temp/start > 1e-6) && (p <= 10000)) {
    p++;
    temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
    start += temp;
  }
  return (4 * M_PI * M_PI * start);
}



// [[Rcpp::export]]
double d_const_vmsin_lambda(double k1, double k2, double lambda)
{
  double result = 0;
  if(lambda != 0){
    double start = 0, rat = pow(lambda, 2) / (4*k1*k2);
    int p = 1;
    double temp = Rf_choose(2*p, p) * p * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
    start += temp;
    while((temp/start > 1e-6) && (p <= 10000)) {
      p++;
      temp = Rf_choose(2*p, p) * p * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
      start += temp;
    }
    result = (8/lambda) * M_PI * M_PI * start;
  }
  return (result);
}

// [[Rcpp::export]]
double d_const_vmsin_k1(double k1, double k2, double lambda)
{
  double start = BESSI(1,k1) * BESSI(0,k2) * 2;
  if(lambda != 0){
    double rat = pow(lambda, 2) / (4*k1*k2);
    int p = 1;
    double temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k2) * ( BESSI((p+1), k1) + BESSI((p-1), k1) - (2*p/k1)*BESSI(p, k1));
    start += temp;
    while((temp/start > 1e-6) && (p <= 10000)) {
      p++;
      temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k2) * ( BESSI((p+1), k1) + BESSI((p-1), k1) - (2*p/k1)*BESSI(p, k1));
      start += temp;
    }
  }
  return (2 * M_PI * M_PI * start);
}

// [[Rcpp::export]]
double d_const_vmsin_k2(double k1, double k2, double lambda)
{
  return (d_const_vmsin_k1(k2, k1, lambda));
}



// [[Rcpp::export]]
arma::vec d_const_vmsin(arma::vec par)
{
  double k1 = par[0], k2 = par[1], lambda = par[2];
  arma::vec result(3);
  result[0] = d_const_vmsin_k1(k1, k2, lambda);
  result[1] = d_const_vmsin_k2(k1, k2, lambda);
  result[2] = d_const_vmsin_lambda(k1, k2, lambda);
  return (result);
}


// [[Rcpp::export]]
arma::vec log_const_vmsin_all(arma::mat par_mat) {
  int K = par_mat.n_cols;
  arma::vec all_lconsts(K);
  for(int j = 0; j < K; j++)
    all_lconsts[j] = log(const_vmsin(par_mat(0,j), par_mat(1,j), par_mat(2,j)));
  return all_lconsts;
}


// [[Rcpp::export]]
double ldsinnum(double x, double y, arma::vec par)
{
  double k1 = par[0], k2 = par[1], lambda = par[2], mu = par[3], nu = par[4];
  return (k1*cos(x-mu) + k2*cos(y-nu) + lambda*sin(x-mu)*sin(y-nu));
}



// [[Rcpp::export]]
arma::vec dsin_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                           arma::vec mu1, arma::vec mu2)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec l_const_all = log_const_vmsin_all(all_par);
  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x[0], x[1], all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec dsin_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                            double mu1, double mu2)
{
  int n = x.n_rows;

  double l_const = log(const_vmsin(k1, k2, k3));
  arma::vec par(5);
  par[0] = k1; par[1] = k2; par[2] = k3; par[3] = mu1; par[4] = mu2;

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), par);
  }

  return arma::exp(ld_num - l_const);
}

// [[Rcpp::export]]
arma::vec dsin_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                             arma::vec mu1, arma::vec mu2)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec l_const_all = log_const_vmsin_all(all_par);

  arma::vec ld_num(n);
   for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}


// generates a single observation
arma::rowvec rsin_single_onepar(double k1, double k2, double k3,
                                double mu1, double mu2, double I_upper_bd) {
  arma::vec par(5);
  par << k1 << k2 << k3 << mu1 << mu2 << arma::endr;

  double x, y, U1, sin_x_minus_mu1, a;

  // First generate x from the marginal distribution given in Eqn. 2.2 in Singh et al 2002

  int accpt_x = 0, rej = 0;
  while (accpt_x < 1) {
    x = runivm_single_onepar(k1, mu1);
    sin_x_minus_mu1 = sin(x - mu1);
    a = sqrt(k2 * k2 + k3 * k3 * sin_x_minus_mu1 * sin_x_minus_mu1);
    U1 = R::unif_rand();
    // I_upper_bd is numerically calculated in R
    if (U1 <= BESSI0(a) / I_upper_bd)
      accpt_x++;
    else
      rej++;
    if(rej == 50) I_upper_bd = BESSI0(sqrt(k2 * k2 + k3 * k3));
    // if there are lots of rejections, use the conservative upper bound.
    if(rej == 5e4) warning("5e4 rejetions.. continuing..");
  }

  double beta = atan(k3 / k2 * sin_x_minus_mu1);

  // Now generate y from vM(a, mu2 + beta) (Eqn. 2.4, Singh et al. 2002)
  y = runivm_single_onepar(a, (mu2 + beta));

  arma::rowvec final_sample(2);
  final_sample[0] = x;
  final_sample[1] = y;

  return final_sample;
}

//  [[Rcpp::export]]
arma::mat rsin_onepar(int n, double k1, double k2, double k3,
                                 double mu1, double mu2, double I_upper_bd) {
  if(n == 1){
    return(rsin_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd));
  } else {
    arma::mat all_sample(n, 2);
    for(int i = 0; i < n; i++)
      all_sample.row(i) = rsin_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd);
    return(all_sample);
  }
}


//  [[Rcpp::export]]
arma::mat rsin_manypar(arma::vec k1, arma::vec k2, arma::vec k3,
                       arma::vec mu1, arma::vec mu2, arma::vec I_upper_bd) {
 int n = k1.size();
    arma::mat all_sample(n, 2);
    for(int i = 0; i < n; i++)
      all_sample.row(i) = rsin_single_onepar(k1[i], k2[i], k3[i], mu1[i], mu2[i], I_upper_bd[i]);
    return(all_sample);
}



// [[Rcpp::export]]
arma::mat mem_p_sin(arma::mat data, arma::mat par, arma::vec pi,
                    arma::vec log_c_von, int ncores = 1)
{
  int n = data.n_rows, K = par.n_cols, i, j;
  double row_total;
  arma::mat den(n, K);
#pragma omp parallel for private(j, row_total) shared(i) num_threads(ncores)
  for(i = 0; i < n; i++){
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = pi[j]*exp(ldsinnum(data(i, 0), data(i, 1), par.col(j)) - log_c_von[j]);;
      row_total += den(i, j);
    }
    row_total = maxi(row_total, 1e-50);
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);
}


// [[Rcpp::export]]
double llik_vmsin_full(arma::mat data, arma::mat par, arma::vec pi, arma::vec log_c, int ncores = 1)
{

  int n = data.n_rows, K = pi.size(), i ,j;
  long double temp, log_sum = 0.0;
  arma::vec log_pi = log(pi);

  if(K > 1) {
#pragma omp parallel for reduction(+:log_sum)  private(j, temp) num_threads(ncores)
    for(i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(ldsinnum(data(i,0), data(i,1), par.col(j)) - log_c[j] + log_pi[j] );
      log_sum += log(temp);
    }
  }  else {
#pragma omp parallel for reduction(+:log_sum) shared(i)  num_threads(ncores)
    for(i = 0; i < n; i++)
      log_sum += ldsinnum(data(i,0), data(i,1), par);
    log_sum -= n*log_c[0];
  }
  return(log_sum);

}


arma::vec grad_log_vmsin_one_comp_i(double x, double y, arma::vec par,
                                    double c_vmsin, arma::vec del_const_vmsin) {
  double k1 = par[0], k2 = par[1], lambda = par[2], mu1 = par[3], mu2 = par[4];
  arma::vec all_entries(5);
  all_entries[0] = (cos(x - mu1) - del_const_vmsin[0]/c_vmsin);
  all_entries[1] = (cos(y - mu2) - del_const_vmsin[1]/c_vmsin);
  all_entries[2] = (sin(x - mu1)*sin(y - mu2) - del_const_vmsin[2]/c_vmsin);
  all_entries[3] = (k1*sin(x - mu1) - lambda*cos(x - mu1)*sin(y - mu2));
  all_entries[4] = (k2*sin(y - mu2) - lambda*sin(x - mu1)*cos(y - mu2));

  return all_entries;
}

// [[Rcpp::export]]
arma::mat grad_vmsin_all_comp(arma::mat data, arma::mat par_mat, arma::vec pi, int ncores = 1)
{
  int n = data.n_rows, K = pi.size(), i ,j;
  arma::vec l_c_vmsin = log_const_vmsin_all(par_mat);
  arma::vec c_vmsin = arma::exp(l_c_vmsin), log_pi = log(pi);

  stdvec_vec del_const_vmsin_all(K);
  stdvec_vec par_mat_cols(K);

  for(j = 0; j < K; j++) {
    del_const_vmsin_all[j] = d_const_vmsin(par_mat.col(j));
    par_mat_cols[j] = par_mat.col(j);
  }

  stdvec_mat grad_sum_part(ncores);
  for(int g = 0; g < ncores; g++)
    grad_sum_part[g] = arma::zeros(5, K);

  double denom, temp;

#pragma omp parallel for   private(j, denom, temp) num_threads(ncores)
  for(i = 0; i < n; i++) {
    arma::mat grad_temp(5, K);
    int g = omp_get_thread_num();
    // int g = 0;
    denom = 0;
    for(j = 0; j < K; j++) {
      temp = exp(ldsinnum(data(i, 0), data(i, 1), par_mat_cols[j]) - l_c_vmsin[j] + log_pi[j]);
      denom += temp;
      grad_temp.col(j) = temp*grad_log_vmsin_one_comp_i(data(i, 0), data(i, 1), par_mat_cols[j], c_vmsin[j], del_const_vmsin_all[j]);
    }
    grad_sum_part[g] += grad_temp/denom;
  }

  arma::mat grad_sum = arma::zeros(5, K);
  for(int g = 0; g < ncores; g++)
    grad_sum += grad_sum_part[g];

  return (grad_sum);
}




// [[Rcpp::export]]
double vmsinmix(double x, double y, arma::mat par, arma::vec pi, arma::vec log_c_von)
{
  double res = 0;
  int K = par.n_cols;

  for(int j = 0; j < K; j++)
    res += exp(ldsinnum(x, y, par.col(j)) - log_c_von[j]) * pi[j];
  return(res);
} //likelihood contribution (density) of a single point in a mixture model

// [[Rcpp::export]]
arma::vec vmsinmix_manyx(arma::mat x, arma::mat par, arma::vec pi, arma::vec log_c_von)
{
  int n = x.n_rows;
  arma::vec result(n);
  for(int i = 0; i < n; i++)
    result[i] = vmsinmix(x(i,0), x(i,1), par, pi, log_c_von);
  return(result);
}
