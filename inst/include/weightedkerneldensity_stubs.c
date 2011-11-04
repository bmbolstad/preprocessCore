#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifndef WEIGHTEDKERNELDENSITY_STUBS_H
#define WEIGHTEDKERNELDENSITY_STUBS_H 1


void KernelDensity(double *x, int *nxxx, double *weights, double *output, double *output_x, int *nout, int *kernel_fn, int *bandwidth_fn, double *bandwidth_adj){
   
  static void(*fun)(double*, int*,  double*, double*, double*, int *, int *, int *, double *) = NULL;
  
  if (fun == NULL)
    fun =  (void(*)(double*, int*,  double*, double*, double*, int *, int *, int *, double *))R_GetCCallable("preprocessCore","KernelDensity");
  
  fun(x, nxxx, weights, output, output_x, nout, kernel_fn, bandwidth_fn, bandwidth_adj);
  
  return;


}



#endif
