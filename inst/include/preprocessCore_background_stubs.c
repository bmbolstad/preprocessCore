#ifndef PREPROCESSCORE_BACKGROUND_STUBS_H
#define PREPROCESSCORE_BACKGROUND_STUBS_H 1


void rma_bg_parameters(double *PM,double *param, int rows, int cols, int column){

  static void(*fun)(double *, double *, int, int, int) = NULL;

  if (fun == NULL)
    fun = (int(*)(double *, double *, int, int, int))R_GetCCallable("preprocessCore","rma_bg_parameters");

  fun(x, param, rows, cols, column);
  return;
}


void rma_bg_adjust(double *PM,double *param, int rows, int cols, int column){

  static void(*fun)(double *, double *, int, int, int) = NULL;

  if (fun == NULL)
    fun = (int(*)(double *, double *, int, int, int))R_GetCCallable("preprocessCore","rma_bg_adjust");

  fun(x, param, rows, cols, column);
  return;
}


void rma_bg_correct(double *PM, int rows, int cols){

  static void(*fun)(double *, int, int) = NULL;

  if (fun == NULL)
    fun = (int(*)(double *, int, int))R_GetCCallable("preprocessCore","rma_bg_correct");

  fun(x, rows, cols);
  return;
}











#endif
