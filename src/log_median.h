#ifndef LOG_MEDIAN_H
#define LOG_MEDIAN_H 1



void LogMedian(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void LogMedian_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

void logmedian(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void logmedian_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);

#endif
