#ifndef BIWEIGHT_H
#define BIWEIGHT_H 1

void tukeybiweight(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void tukeybiweight_no_log(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void TukeyBiweight(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void TukeyBiweight_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);
void TukeyBiweight_no_log_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

double Tukey_Biweight(double *x, size_t length);

#endif
