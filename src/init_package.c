/*****************************************************
 **
 ** file: init_package.c
 **
 ** Copyright (C) 2007    B. M. Bolstad
 **
 ** aim: Register c code routines so that they can be called in other packages.
 **
 ** History
 ** May 20, 2007 - Initial version
 **
 *****************************************************/

#include "qnorm.h"
#include "medianpolish.h"
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RConverters.h>

#if _MSC_VER >= 1000
__declspec(dllexport)
#endif


static const R_CallMethodDef callMethods[]  = {
  {"R_qnorm_c",(DL_FUNC)&R_qnorm_c,2},
    {"R_qnorm_robust_weights", (DL_FUNC)&R_qnorm_robust_weights, 3},
    {"R_qnorm_robust_c",(DL_FUNC)&R_qnorm_robust_c,6},
    {"R_qnorm_determine_target",(DL_FUNC)&R_qnorm_determine_target,2},
    {"R_qnorm_using_target",(DL_FUNC)&R_qnorm_using_target,3 },
    {"R_qnorm_within_blocks",(DL_FUNC)&R_qnorm_within_blocks, 3},
    {NULL, NULL, 0}
  };

void R_init_preprocessCore(DllInfo *info){

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);


  /* The normalization routines */
  R_RegisterCCallable("preprocessCore", "qnorm_c", (DL_FUNC)&qnorm_c);
  R_RegisterCCallable("preprocessCore", "qnorm_robust_c", (DL_FUNC)&qnorm_robust_c);
  R_RegisterCCallable("preprocessCore", "qnorm_c_using_target", (DL_FUNC)&qnorm_c_using_target);
  R_RegisterCCallable("preprocessCore", "qnorm_c_determine_target", (DL_FUNC)&qnorm_c_determine_target);
  R_RegisterCCallable("preprocessCore", "qnorm_c_within_blocks", (DL_FUNC)&qnorm_c_within_blocks);

  /* The summarization routines */

  R_RegisterCCallable("preprocessCore", "median_polish_fit_no_copy", (DL_FUNC)&median_polish_fit_no_copy);
  R_RegisterCCallable("preprocessCore", "median_polish_no_copy", (DL_FUNC)&median_polish_no_copy);
  R_RegisterCCallable("preprocessCore", "median_polish_log2_no_copy", (DL_FUNC)&median_polish_log2_no_copy);
  R_RegisterCCallable("preprocessCore", "median_polish_log2", (DL_FUNC)&median_polish_log2);
  R_RegisterCCallable("preprocessCore", "median_polish", (DL_FUNC)&median_polish);
  R_RegisterCCallable("preprocessCore", "MedianPolish", (DL_FUNC)&MedianPolish);
  


}
