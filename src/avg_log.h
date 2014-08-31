#ifndef AVG_LOG_H
#define AVG_LOG_H

void AverageLog(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE);
void AverageLog_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes);

void averagelog_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE);
void averagelog(double *data, size_t rows, size_t cols, double *results, double *resultsSE);

#endif
