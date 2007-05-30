/*****************************************************
 **
 ** file: init_package.c
 **
 ** Copyright (C) 2007    B. M. Bolstad
 **
 ** aim: Register c code routines so that they can be called in other packages.
 **"
 ** History
 ** May 20, 2007 - Initial version
 ** May 24-27, 2007 - add in additional registered functions
 **
 *****************************************************/

#include "qnorm.h"
#include "medianpolish.h"
#include "log_avg.h"
#include "avg_log.h"
#include "biweight.h"
#include "lm.h"
#include "rlm.h"
#include "rlm_se.h"

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
  

  R_RegisterCCallable("preprocessCore","AverageLog", (DL_FUNC)&AverageLog);
  R_RegisterCCallable("preprocessCore","averagelog_no_copy", (DL_FUNC)&averagelog_no_copy);
  R_RegisterCCallable("preprocessCore","averagelog", (DL_FUNC)&averagelog);
  R_RegisterCCallable("preprocessCore","AverageLog_noSE", (DL_FUNC)&AverageLog_noSE);

  R_RegisterCCallable("preprocessCore","logaverage", (DL_FUNC)&logaverage);
  R_RegisterCCallable("preprocessCore","LogAverage", (DL_FUNC)&LogAverage);

  R_RegisterCCallable("preprocessCore","tukeybiweight", (DL_FUNC)&tukeybiweight);
  R_RegisterCCallable("preprocessCore","TukeyBiweight", (DL_FUNC)&TukeyBiweight);
  R_RegisterCCallable("preprocessCore","Tukey_Biweight", (DL_FUNC)&Tukey_Biweight);


  R_RegisterCCallable("preprocessCore","lm_wfit", (DL_FUNC)&lm_wfit);
  
  
  R_RegisterCCallable("preprocessCore","rlm_fit", (DL_FUNC)&rlm_fit);
  R_RegisterCCallable("preprocessCore","rlm_wfit", (DL_FUNC)&rlm_wfit);

  R_RegisterCCallable("preprocessCore","rlm_compute_se", (DL_FUNC)&rlm_compute_se);

  
  R_RegisterCCallable("preprocessCore","rlm_fit_anova", (DL_FUNC)&rlm_fit_anova);
  R_RegisterCCallable("preprocessCore","rlm_wfit_anova", (DL_FUNC)&rlm_wfit_anova);
  R_RegisterCCallable("preprocessCore","rlm_compute_se_anova", (DL_FUNC)&rlm_compute_se_anova);
   


}
