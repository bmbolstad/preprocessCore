#ifndef WEIGHTEDKERNELDENSITY_H
#define WEIGHTEDKERNELDENSITY_H

void KernelDensity(double *x, size_t nxxx, double *weights, double *output, double *output_x, size_t nout, int kernel_fn, int bandwidth_fn, double bandwidth_adj);
void KernelDensity_lowmem(double *x, size_t nxxx, double *output, double *output_x, size_t nout);

#endif
