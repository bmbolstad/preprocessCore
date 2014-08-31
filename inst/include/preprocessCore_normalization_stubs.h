int qnorm_c(double *data, int *rows, int *cols);
int qnorm_robust_c(double *data,double *weights, int *rows, int *cols, int *use_median, int *use_log2, int *weight_scheme);
int qnorm_c_using_target(double *data, int *rows, int *cols, double *target, int *targetrows);
int qnorm_c_determine_target(double *data, int *rows, int *cols, double *target, int *targetrows);
int qnorm_c_within_blocks(double *x, int *rows, int *cols, int *blocks);
