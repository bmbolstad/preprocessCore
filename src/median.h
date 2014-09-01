#ifndef MEDIAN_H
#define MEDIAN_H 1



void ColMedian(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void ColMedian_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

void colmedian(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void colmedian_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);




#endif
