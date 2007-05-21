#ifndef QNORM_H
#define QNORM_H 1

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
 


int qnorm_c(double *data, int *rows, int *cols);
int qnorm_robust_c(double *data,double *weights, int *rows, int *cols, int *use_median, int *use_log2, int *weight_scheme);

SEXP R_qnorm_robust_weights(SEXP X, SEXP remove_extreme, SEXP n_remove);

#endif


int qnorm_d(double *data, int *rows, int *cols);
