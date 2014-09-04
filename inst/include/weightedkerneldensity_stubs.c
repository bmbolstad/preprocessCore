#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifndef WEIGHTEDKERNELDENSITY_STUBS_H
#define WEIGHTEDKERNELDENSITY_STUBS_H 1

#include <stddef.h>

void KernelDensity(double *x, size_t nxxx, double *weights, double *output, double *output_x, size_t nout, int kernel_fn, int bandwidth_fn, double bandwidth_adj){
   
  static void(*fun)(double*, size_t,  double*, double*, double*, size_t, int, int, double) = NULL;
  
  if (fun == NULL)
    fun =  (void(*)(double*, size_t,  double*, double*, double*, size_t, int, int, double))R_GetCCallable("preprocessCore","KernelDensity");
  
  fun(x, nxxx, weights, output, output_x, nout, kernel_fn, bandwidth_fn, bandwidth_adj);
  
  return;


}



#endif
