#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BAMBI_BESSI0_C(SEXP);
extern SEXP _BAMBI_cID(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_const_univm(SEXP);
extern SEXP _BAMBI_const_uniwnorm(SEXP);
extern SEXP _BAMBI_const_vmcos(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_const_vmcos_anltc(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_const_vmcos_mc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_const_vmsin(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_const_wnorm2(SEXP);
extern SEXP _BAMBI_d_const_vmcos(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmcos_anltc(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmcos_k1_anltc(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmcos_mc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmsin(SEXP);
extern SEXP _BAMBI_d_const_vmsin_k1(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmsin_k2(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_d_const_vmsin_lambda(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dcos_manyx_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dcos_manyx_onepar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dcos_onex_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dsin_manyx_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dsin_manyx_onepar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dsin_onex_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dunivm_manyx_manypar(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dunivm_manyx_onepar(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dunivm_onex_manypar(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_duniwnorm_manyx_manypar(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_duniwnorm_manyx_onepar(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_duniwnorm_onex_manypar(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dwnorm2_manyx_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dwnorm2_manyx_onepar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_dwnorm2_onex_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_den_wnorm2_one_comp_i_unadj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_log_vmcos_one_comp_i(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_univm_all_comp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_uniwnorm_all_comp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_vmcos_all_comp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_vmsin_all_comp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_grad_wnorm2_all_comp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_l_const_uniwnorm(SEXP);
extern SEXP _BAMBI_l_const_wnorm2(SEXP);
extern SEXP _BAMBI_ldcosnum(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_ldsinnum(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_ldunivmnum(SEXP, SEXP);
extern SEXP _BAMBI_lduniwnormnum(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_ldwnorm2_num(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_llik_univm_full(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_llik_uniwnorm_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_llik_vmcos_full(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_llik_vmsin_full(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_llik_wnorm2_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_log_const_univm_all(SEXP);
extern SEXP _BAMBI_log_const_uniwnorm_all(SEXP);
extern SEXP _BAMBI_log_const_vmcos_all(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_log_const_vmsin_all(SEXP);
extern SEXP _BAMBI_log_const_wnorm2_all(SEXP);
extern SEXP _BAMBI_mem_p_cos(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_mem_p_sin(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_mem_p_univm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_mem_p_uniwnorm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_mem_p_wnorm2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_par_mat_permute(SEXP, SEXP);
extern SEXP _BAMBI_rcos_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_rcos_onepar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_rowVars(SEXP);
extern SEXP _BAMBI_rsin_manypar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_rsin_onepar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_runivm_manypar(SEXP, SEXP);
extern SEXP _BAMBI_runivm_onepar(SEXP, SEXP, SEXP);
extern SEXP _BAMBI_univmmix(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_univmmix_manyx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_uniwnormmix(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_uniwnormmix_manyx(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_vmcosmix(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_vmcosmix_manyx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_vmsinmix(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_vmsinmix_manyx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_wnorm2mix(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BAMBI_wnorm2mix_manyx(SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
  {"_BAMBI_BESSI0_C",                         (DL_FUNC) &_BAMBI_BESSI0_C,                         1},
  {"_BAMBI_cID",                              (DL_FUNC) &_BAMBI_cID,                              3},
  {"_BAMBI_const_univm",                      (DL_FUNC) &_BAMBI_const_univm,                      1},
  {"_BAMBI_const_uniwnorm",                   (DL_FUNC) &_BAMBI_const_uniwnorm,                   1},
  {"_BAMBI_const_vmcos",                      (DL_FUNC) &_BAMBI_const_vmcos,                      5},
  {"_BAMBI_const_vmcos_anltc",                (DL_FUNC) &_BAMBI_const_vmcos_anltc,                3},
  {"_BAMBI_const_vmcos_mc",                   (DL_FUNC) &_BAMBI_const_vmcos_mc,                   5},
  {"_BAMBI_const_vmsin",                      (DL_FUNC) &_BAMBI_const_vmsin,                      3},
  {"_BAMBI_const_wnorm2",                     (DL_FUNC) &_BAMBI_const_wnorm2,                     1},
  {"_BAMBI_d_const_vmcos",                    (DL_FUNC) &_BAMBI_d_const_vmcos,                    3},
  {"_BAMBI_d_const_vmcos_anltc",              (DL_FUNC) &_BAMBI_d_const_vmcos_anltc,              3},
  {"_BAMBI_d_const_vmcos_k1_anltc",           (DL_FUNC) &_BAMBI_d_const_vmcos_k1_anltc,           3},
  {"_BAMBI_d_const_vmcos_mc",                 (DL_FUNC) &_BAMBI_d_const_vmcos_mc,                 5},
  {"_BAMBI_d_const_vmsin",                    (DL_FUNC) &_BAMBI_d_const_vmsin,                    1},
  {"_BAMBI_d_const_vmsin_k1",                 (DL_FUNC) &_BAMBI_d_const_vmsin_k1,                 3},
  {"_BAMBI_d_const_vmsin_k2",                 (DL_FUNC) &_BAMBI_d_const_vmsin_k2,                 3},
  {"_BAMBI_d_const_vmsin_lambda",             (DL_FUNC) &_BAMBI_d_const_vmsin_lambda,             3},
  {"_BAMBI_dcos_manyx_manypar",               (DL_FUNC) &_BAMBI_dcos_manyx_manypar,               7},
  {"_BAMBI_dcos_manyx_onepar",                (DL_FUNC) &_BAMBI_dcos_manyx_onepar,                7},
  {"_BAMBI_dcos_onex_manypar",                (DL_FUNC) &_BAMBI_dcos_onex_manypar,                7},
  {"_BAMBI_dsin_manyx_manypar",               (DL_FUNC) &_BAMBI_dsin_manyx_manypar,               6},
  {"_BAMBI_dsin_manyx_onepar",                (DL_FUNC) &_BAMBI_dsin_manyx_onepar,                6},
  {"_BAMBI_dsin_onex_manypar",                (DL_FUNC) &_BAMBI_dsin_onex_manypar,                6},
  {"_BAMBI_dunivm_manyx_manypar",             (DL_FUNC) &_BAMBI_dunivm_manyx_manypar,             3},
  {"_BAMBI_dunivm_manyx_onepar",              (DL_FUNC) &_BAMBI_dunivm_manyx_onepar,              3},
  {"_BAMBI_dunivm_onex_manypar",              (DL_FUNC) &_BAMBI_dunivm_onex_manypar,              3},
  {"_BAMBI_duniwnorm_manyx_manypar",          (DL_FUNC) &_BAMBI_duniwnorm_manyx_manypar,          4},
  {"_BAMBI_duniwnorm_manyx_onepar",           (DL_FUNC) &_BAMBI_duniwnorm_manyx_onepar,           4},
  {"_BAMBI_duniwnorm_onex_manypar",           (DL_FUNC) &_BAMBI_duniwnorm_onex_manypar,           4},
  {"_BAMBI_dwnorm2_manyx_manypar",            (DL_FUNC) &_BAMBI_dwnorm2_manyx_manypar,            7},
  {"_BAMBI_dwnorm2_manyx_onepar",             (DL_FUNC) &_BAMBI_dwnorm2_manyx_onepar,             7},
  {"_BAMBI_dwnorm2_onex_manypar",             (DL_FUNC) &_BAMBI_dwnorm2_onex_manypar,             7},
  {"_BAMBI_grad_den_wnorm2_one_comp_i_unadj", (DL_FUNC) &_BAMBI_grad_den_wnorm2_one_comp_i_unadj, 6},
  {"_BAMBI_grad_log_vmcos_one_comp_i",        (DL_FUNC) &_BAMBI_grad_log_vmcos_one_comp_i,        5},
  {"_BAMBI_grad_univm_all_comp",              (DL_FUNC) &_BAMBI_grad_univm_all_comp,              4},
  {"_BAMBI_grad_uniwnorm_all_comp",           (DL_FUNC) &_BAMBI_grad_uniwnorm_all_comp,           5},
  {"_BAMBI_grad_vmcos_all_comp",              (DL_FUNC) &_BAMBI_grad_vmcos_all_comp,              5},
  {"_BAMBI_grad_vmsin_all_comp",              (DL_FUNC) &_BAMBI_grad_vmsin_all_comp,              4},
  {"_BAMBI_grad_wnorm2_all_comp",             (DL_FUNC) &_BAMBI_grad_wnorm2_all_comp,             5},
  {"_BAMBI_l_const_uniwnorm",                 (DL_FUNC) &_BAMBI_l_const_uniwnorm,                 1},
  {"_BAMBI_l_const_wnorm2",                   (DL_FUNC) &_BAMBI_l_const_wnorm2,                   1},
  {"_BAMBI_ldcosnum",                         (DL_FUNC) &_BAMBI_ldcosnum,                         3},
  {"_BAMBI_ldsinnum",                         (DL_FUNC) &_BAMBI_ldsinnum,                         3},
  {"_BAMBI_ldunivmnum",                       (DL_FUNC) &_BAMBI_ldunivmnum,                       2},
  {"_BAMBI_lduniwnormnum",                    (DL_FUNC) &_BAMBI_lduniwnormnum,                    3},
  {"_BAMBI_ldwnorm2_num",                     (DL_FUNC) &_BAMBI_ldwnorm2_num,                     3},
  {"_BAMBI_llik_univm_full",                  (DL_FUNC) &_BAMBI_llik_univm_full,                  5},
  {"_BAMBI_llik_uniwnorm_full",               (DL_FUNC) &_BAMBI_llik_uniwnorm_full,               6},
  {"_BAMBI_llik_vmcos_full",                  (DL_FUNC) &_BAMBI_llik_vmcos_full,                  5},
  {"_BAMBI_llik_vmsin_full",                  (DL_FUNC) &_BAMBI_llik_vmsin_full,                  5},
  {"_BAMBI_llik_wnorm2_full",                 (DL_FUNC) &_BAMBI_llik_wnorm2_full,                 6},
  {"_BAMBI_log_const_univm_all",              (DL_FUNC) &_BAMBI_log_const_univm_all,              1},
  {"_BAMBI_log_const_uniwnorm_all",           (DL_FUNC) &_BAMBI_log_const_uniwnorm_all,           1},
  {"_BAMBI_log_const_vmcos_all",              (DL_FUNC) &_BAMBI_log_const_vmcos_all,              3},
  {"_BAMBI_log_const_vmsin_all",              (DL_FUNC) &_BAMBI_log_const_vmsin_all,              1},
  {"_BAMBI_log_const_wnorm2_all",             (DL_FUNC) &_BAMBI_log_const_wnorm2_all,             1},
  {"_BAMBI_mem_p_cos",                        (DL_FUNC) &_BAMBI_mem_p_cos,                        5},
  {"_BAMBI_mem_p_sin",                        (DL_FUNC) &_BAMBI_mem_p_sin,                        5},
  {"_BAMBI_mem_p_univm",                      (DL_FUNC) &_BAMBI_mem_p_univm,                      5},
  {"_BAMBI_mem_p_uniwnorm",                   (DL_FUNC) &_BAMBI_mem_p_uniwnorm,                   6},
  {"_BAMBI_mem_p_wnorm2",                     (DL_FUNC) &_BAMBI_mem_p_wnorm2,                     6},
  {"_BAMBI_par_mat_permute",                  (DL_FUNC) &_BAMBI_par_mat_permute,                  2},
  {"_BAMBI_rcos_manypar",                     (DL_FUNC) &_BAMBI_rcos_manypar,                     6},
  {"_BAMBI_rcos_onepar",                      (DL_FUNC) &_BAMBI_rcos_onepar,                      7},
  {"_BAMBI_rowVars",                          (DL_FUNC) &_BAMBI_rowVars,                          1},
  {"_BAMBI_rsin_manypar",                     (DL_FUNC) &_BAMBI_rsin_manypar,                     6},
  {"_BAMBI_rsin_onepar",                      (DL_FUNC) &_BAMBI_rsin_onepar,                      7},
  {"_BAMBI_runivm_manypar",                   (DL_FUNC) &_BAMBI_runivm_manypar,                   2},
  {"_BAMBI_runivm_onepar",                    (DL_FUNC) &_BAMBI_runivm_onepar,                    3},
  {"_BAMBI_univmmix",                         (DL_FUNC) &_BAMBI_univmmix,                         4},
  {"_BAMBI_univmmix_manyx",                   (DL_FUNC) &_BAMBI_univmmix_manyx,                   4},
  {"_BAMBI_uniwnormmix",                      (DL_FUNC) &_BAMBI_uniwnormmix,                      5},
  {"_BAMBI_uniwnormmix_manyx",                (DL_FUNC) &_BAMBI_uniwnormmix_manyx,                5},
  {"_BAMBI_vmcosmix",                         (DL_FUNC) &_BAMBI_vmcosmix,                         5},
  {"_BAMBI_vmcosmix_manyx",                   (DL_FUNC) &_BAMBI_vmcosmix_manyx,                   4},
  {"_BAMBI_vmsinmix",                         (DL_FUNC) &_BAMBI_vmsinmix,                         5},
  {"_BAMBI_vmsinmix_manyx",                   (DL_FUNC) &_BAMBI_vmsinmix_manyx,                   4},
  {"_BAMBI_wnorm2mix",                        (DL_FUNC) &_BAMBI_wnorm2mix,                        5},
  {"_BAMBI_wnorm2mix_manyx",                  (DL_FUNC) &_BAMBI_wnorm2mix_manyx,                  5},
  {NULL, NULL, 0}
};

void R_init_BAMBI(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
