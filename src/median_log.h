#ifndef MEDIAN_LOG_H
#define MEDIAN_LOG_H 1



void MedianLog(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void MedianLog_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

void medianlog(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void medianlog_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);




#endif
