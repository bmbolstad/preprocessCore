#ifndef MEDIANPOLISH_H
#define MEDIANPOLISH_H 1

void median_polish_fit_no_copy(double *data, size_t rows, size_t cols, double *r, double *c, double *t);
void median_polish_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void median_polish_log2_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void median_polish_log2(double *data, size_t rows, size_t cols, double *results, double *resultsSE, double *residuals);
void median_polish(double *data, size_t rows, size_t cols, double *results, double *resultsSE, double *residuals);
void MedianPolish(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void MedianPolish_no_log(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);


#endif
