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
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>



void R_init_preprocessCore(){


  R_RegisterCCallable("preprocessCore", "qnorm_c", (DL_FUNC)&qnorm_c);


}
