#ifndef AVG_H
#define AVG_H

void ColAverage(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void ColAverage_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

void colaverage_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void colaverage(double *data, size_t rows, size_t cols, double *results, double *resultsSE);

#endif
