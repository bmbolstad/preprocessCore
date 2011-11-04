#ifndef WEIGHTEDKERNELDENSITY_H
#define WEIGHTEDKERNELDENSITY_H

void KernelDensity(double *x, int *nxxx, double *weights, double *output, double *output_x, int *nout, int *kernel_fn, int *bandwidth_fn, double *bandwidth_adj);
void KernelDensity_lowmem(double *x, int *nxxx, double *output, double *output_x, int *nout);

#endif
