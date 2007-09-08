#ifndef BIWEIGHT_H
#define BIWEIGHT_H 1

void tukeybiweight(double *data, int rows, int cols, double *results, double *resultsSE);
void TukeyBiweight(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE);
double Tukey_Biweight(double *x, int length);

#endif
