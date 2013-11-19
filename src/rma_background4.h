#ifndef RMA_BACKGROUND4_H
#define RMA_BACKGROUND4_H

void rma_bg_parameters(double *PM,double *param, size_t rows, size_t cols, size_t column);
void rma_bg_adjust(double *PM, double *param, size_t rows, size_t cols, size_t column);
void rma_bg_correct(double *PM, size_t rows, size_t cols);

SEXP R_rma_bg_correct(SEXP PMmat,SEXP copy);

#endif
