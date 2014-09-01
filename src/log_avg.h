#ifndef LOG_AVG_H
#define LOG_AVG_H 1

void logaverage(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void LogAverage(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void LogAverage_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

#endif
