/*********************************************************************
 **
 ** file: R_subrcModel_interface.c
 **
 ** Aim: Code which provides .Call() interfaces to the subset of rows in 
 ** a matrix rcModel fitting
 **
 ** Copyright (C) 2012 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bmb@bmbolstad.com>
 ** 
 ** created on: Mar 7, 2012
 **
 ** History
 ** Mar 7, 2012 - Initial version
 **
 *********************************************************************/


#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "rlm.h"
#include "rlm_se.h"
#include "psi_fns.h"
#include "medianpolish.h"






#ifdef USE_PTHREADS
#include <pthread.h>
#include <limits.h>
#include <unistd.h>
#define THREADS_ENV_VAR "R_THREADS"
pthread_mutex_t mutex_R;
struct loop_data{
  double *matrix;
  SEXP *R_return_value;
  SEXP *R_rowIndexList;
  SEXP *PsiCode;
  SEXP *PsiK;
  SEXP *Scales; 
  int rows;
  int cols;
  int length_rowIndexList;
  int start_row;
  int end_row;
};
#endif




#ifdef  USE_PTHREADS
static void *sub_rcModelSummarize_medianpolish_group(void *data){
  
  struct loop_data *args = (struct loop_data *) data;
  int *cur_rows;
  double *buffer, *buffer2;
  int i, j, k;
  int ncur_rows;

  
  SEXP R_return_value_cur;
  SEXP R_weights;
  SEXP R_residuals;
  SEXP R_beta;
  SEXP R_SE;
  
  SEXP R_return_value_names;

  double *beta;
  double *residuals;
  double *weights;
  double *se;

  double intercept;

  int cols = args->cols;
 
  for (j = args->start_row; j <= args->end_row;  j++){    
    ncur_rows = LENGTH(VECTOR_ELT(*(args->R_rowIndexList),j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(*(args->R_rowIndexList),j));

    PROTECT(R_return_value_cur = allocVector(VECSXP,4));
    PROTECT(R_beta = allocVector(REALSXP, ncur_rows + cols));
    /* PROTECT(R_weights = allocMatrix(REALSXP,ncur_rows,cols));*/
    PROTECT(R_residuals = allocMatrix(REALSXP,ncur_rows,cols));
    /*  PROTECT(R_SE = allocVector(REALSXP,ncur_rows+cols)); */

    R_weights = R_NilValue;
    R_SE = R_NilValue;


    SET_VECTOR_ELT(R_return_value_cur,0,R_beta);
    SET_VECTOR_ELT(R_return_value_cur,1,R_weights);
    SET_VECTOR_ELT(R_return_value_cur,2,R_residuals);
    SET_VECTOR_ELT(R_return_value_cur,3,R_SE);

    UNPROTECT(2);

    beta = NUMERIC_POINTER(R_beta);
    residuals = NUMERIC_POINTER(R_residuals);
    /*  weights = NUMERIC_POINTER(R_weights);
        se = NUMERIC_POINTER(R_SE);
    */

    for (k = 0; k < cols; k++){
        for (i =0; i < ncur_rows; i++){
     	    residuals[k*ncur_rows + i] = args->matrix[k*args->rows + cur_rows[i]];  
        }
    } 

    memset(beta, 0, (ncur_rows+cols)*sizeof(double));

    median_polish_fit_no_copy(residuals, ncur_rows, cols, &beta[cols], &beta[0], &intercept);

    for (i=0; i < cols; i++)
        beta[i]+=intercept;

    PROTECT(R_return_value_names= allocVector(STRSXP,4));
    SET_STRING_ELT(R_return_value_names,0,mkChar("Estimates"));
    SET_STRING_ELT(R_return_value_names,1,mkChar("Weights"));
    SET_STRING_ELT(R_return_value_names,2,mkChar("Residuals"));
    SET_STRING_ELT(R_return_value_names,3,mkChar("StdErrors"));
    setAttrib(R_return_value_cur, R_NamesSymbol,R_return_value_names);
    UNPROTECT(2);
   
    pthread_mutex_lock (&mutex_R);
    SET_VECTOR_ELT(*(args->R_return_value),j,R_return_value_cur);
    pthread_mutex_unlock (&mutex_R);
  }
}
#endif




SEXP R_sub_rcModelSummarize_medianpolish(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_return_value;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer, *buffer2;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;
#ifdef USE_PTHREADS
  int t, returnCode, chunk_size, num_threads = 1;
  double chunk_size_d, chunk_tot_d;
  char *nthreads;
  pthread_attr_t attr;
  pthread_t *threads;
  struct loop_data *args;
  void *status; 
#ifdef PTHREAD_STACK_MIN
  size_t stacksize = PTHREAD_STACK_MIN + 0x4000;
#else
  size_t stacksize = 0x8000;
#endif
#else

  SEXP R_return_value_cur;

  SEXP R_weights;
  SEXP R_residuals;
  SEXP R_beta;
  SEXP R_SE;
  
  SEXP R_return_value_names;

  double *beta;
  double *residuals;
  double *weights;
  double *se;

  double intercept;
#endif

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_return_value = allocVector(VECSXP,length_rowIndexList));
  
#ifdef  USE_PTHREADS
  nthreads = getenv(THREADS_ENV_VAR);
  if(nthreads != NULL){
    num_threads = atoi(nthreads);
    if(num_threads <= 0){
      error("The number of threads (enviroment variable %s) must be a positive integer, but the specified value was %s", THREADS_ENV_VAR, nthreads);
    }
  }
  threads = (pthread_t *) Calloc(num_threads, pthread_t);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setstacksize (&attr, stacksize);
  
  /* this code works out how many threads to use and allocates ranges of subColumns to each thread */
  /* The aim is to try to be as fair as possible in dividing up the matrix */
  /* A special cases to be aware of: 
     1) Number of subColumns is less than the number of threads
  */
  
  if (num_threads < length_rowIndexList){
    chunk_size = length_rowIndexList/num_threads;
    chunk_size_d = ((double) length_rowIndexList)/((double) num_threads);
  } else {
    chunk_size = 1;
    chunk_size_d = 1;
  }

  if(chunk_size == 0){
    chunk_size = 1;
  }
  args = (struct loop_data *) Calloc((length_rowIndexList < num_threads ? length_rowIndexList : num_threads), struct loop_data);

  args[0].matrix = matrix;
  args[0].R_return_value = &R_return_value;
  args[0].R_rowIndexList = &R_rowIndexList;
  args[0].rows = rows;  
  args[0].cols = cols;
  args[0].length_rowIndexList = length_rowIndexList;

  pthread_mutex_init(&mutex_R, NULL);

  t = 0; /* t = number of actual threads doing work */
  chunk_tot_d = 0;
  for (i=0; floor(chunk_tot_d+0.00001) < length_rowIndexList; i+=chunk_size){
     if(t != 0){
       memcpy(&(args[t]), &(args[0]), sizeof(struct loop_data));
     }

     args[t].start_row = i;     
     /* take care of distribution of the remainder (when #chips%#threads != 0) */
     chunk_tot_d += chunk_size_d;
     // Add 0.00001 in case there was a rounding issue with the division
     if(i+chunk_size < floor(chunk_tot_d+0.00001)){
       args[t].end_row = i+chunk_size;
       i++;
     }
     else{
       args[t].end_row = i+chunk_size-1;
     }
     t++;
  }

  
  for (i =0; i < t; i++){
     returnCode = pthread_create(&threads[i], &attr, sub_rcModelSummarize_medianpolish_group, (void *) &(args[i]));
     if (returnCode){
         error("ERROR; return code from pthread_create() is %d\n", returnCode);
     }
  }
  /* Wait for the other threads */
  for(i = 0; i < t; i++){
      returnCode = pthread_join(threads[i], &status);
      if (returnCode){
         error("ERROR; return code from pthread_join(thread #%d) is %d, exit status for thread was %d\n", 
               i, returnCode, *((int *) status));
      }
  }

  pthread_attr_destroy(&attr);  
  pthread_mutex_destroy(&mutex_R);
  Free(threads);
  Free(args);  
#else     

  for (j =0; j < length_rowIndexList; j++){    

    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));

    PROTECT(R_return_value_cur = allocVector(VECSXP,4));
    PROTECT(R_beta = allocVector(REALSXP, ncur_rows + cols));
    /* PROTECT(R_weights = allocMatrix(REALSXP,ncur_rows,cols));*/
    PROTECT(R_residuals = allocMatrix(REALSXP,ncur_rows,cols));
    /*  PROTECT(R_SE = allocVector(REALSXP,ncur_rows+cols)); */

    R_weights = R_NilValue;
    R_SE = R_NilValue;


    SET_VECTOR_ELT(R_return_value_cur,0,R_beta);
    SET_VECTOR_ELT(R_return_value_cur,1,R_weights);
    SET_VECTOR_ELT(R_return_value_cur,2,R_residuals);
    SET_VECTOR_ELT(R_return_value_cur,3,R_SE);

    UNPROTECT(2);

    beta = NUMERIC_POINTER(R_beta);
    residuals = NUMERIC_POINTER(R_residuals);
    /*  weights = NUMERIC_POINTER(R_weights);
        se = NUMERIC_POINTER(R_SE);
    */

    for (k = 0; k < cols; k++){
        for (i =0; i < ncur_rows; i++){
     	    residuals[k*ncur_rows + i] = matrix[k*rows + cur_rows[i]];  
        }
    } 

    memset(beta, 0, (ncur_rows+cols)*sizeof(double));

    median_polish_fit_no_copy(residuals, ncur_rows, cols, &beta[cols], &beta[0], &intercept);

    for (i=0; i < cols; i++)
        beta[i]+=intercept;

    PROTECT(R_return_value_names= allocVector(STRSXP,4));
    SET_STRING_ELT(R_return_value_names,0,mkChar("Estimates"));
    SET_STRING_ELT(R_return_value_names,1,mkChar("Weights"));
    SET_STRING_ELT(R_return_value_names,2,mkChar("Residuals"));
    SET_STRING_ELT(R_return_value_names,3,mkChar("StdErrors"));
    setAttrib(R_return_value_cur, R_NamesSymbol,R_return_value_names);
    UNPROTECT(2);
    SET_VECTOR_ELT(R_return_value,j,R_return_value_cur);
  }
#endif
  UNPROTECT(1);
  return R_return_value;
}












#ifdef  USE_PTHREADS
static void *sub_rcModelSummarize_plm_group(void *data){
  
  struct loop_data *args = (struct loop_data *) data;
  int *cur_rows;
  double *buffer, *buffer2;
  int i, j, k;
  int ncur_rows;

  
  SEXP R_return_value_cur;
  SEXP R_weights;
  SEXP R_residuals;
  SEXP R_beta;
  SEXP R_SE;
  SEXP R_scale;

  SEXP R_return_value_names;

  double *Ymat;

  double *beta;
  double *residuals;
  double *weights;
  double *se;
  
  double scale=-1.0;
  double *scaleptr;

  double residSE;

  int cols = args->cols;
 
  for (j = args->start_row; j <= args->end_row;  j++){    
    ncur_rows = LENGTH(VECTOR_ELT(*(args->R_rowIndexList),j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(*(args->R_rowIndexList),j));

    PROTECT(R_return_value_cur = allocVector(VECSXP,5));
    PROTECT(R_beta = allocVector(REALSXP, ncur_rows + cols));
    PROTECT(R_weights = allocMatrix(REALSXP,ncur_rows,cols));
    PROTECT(R_residuals = allocMatrix(REALSXP,ncur_rows,cols));
    PROTECT(R_SE = allocVector(REALSXP,ncur_rows+cols)); 
    PROTECT(R_scale = allocVector(REALSXP,1));

    SET_VECTOR_ELT(R_return_value_cur,0,R_beta);
    SET_VECTOR_ELT(R_return_value_cur,1,R_weights);
    SET_VECTOR_ELT(R_return_value_cur,2,R_residuals);
    SET_VECTOR_ELT(R_return_value_cur,3,R_SE);
    SET_VECTOR_ELT(R_return_value_cur,4,R_scale);

    UNPROTECT(5);

    beta = NUMERIC_POINTER(R_beta);
    residuals = NUMERIC_POINTER(R_residuals);
    weights = NUMERIC_POINTER(R_weights);
    se = NUMERIC_POINTER(R_SE);
    

   scaleptr = NUMERIC_POINTER(R_scale);


    if (isNull(*args->Scales)){
      scaleptr[0] = -1.0;
    } else if (length(*args->Scales) != cols) {
      scaleptr[0] = NUMERIC_POINTER(*args->Scales)[0];
    }


    Ymat = Calloc(ncur_rows*cols,double);
    
    
    for (k = 0; k < cols; k++){
        for (i =0; i < ncur_rows; i++){
     	    Ymat[k*ncur_rows + i] = args->matrix[k*args->rows + cur_rows[i]];  
        }
    } 

    rlm_fit_anova_scale(Ymat, ncur_rows, cols, scaleptr, beta, residuals, weights, PsiFunc(asInteger(*args->PsiCode)),asReal(*args->PsiK), 20, 0);
  
    rlm_compute_se_anova(Ymat, ncur_rows, cols, beta, residuals, weights,se, (double *)NULL, &residSE, 4, PsiFunc(asInteger(*args->PsiCode)),asReal(*args->PsiK));


  

     beta[ncur_rows+cols -1] = 0.0;
  

     for (i = cols; i < ncur_rows + cols -1; i++)
        beta[ncur_rows+cols -1]-=beta[i];

     Free(Ymat);
     PROTECT(R_return_value_names= allocVector(STRSXP,5));
     SET_STRING_ELT(R_return_value_names,0,mkChar("Estimates"));
     SET_STRING_ELT(R_return_value_names,1,mkChar("Weights"));
     SET_STRING_ELT(R_return_value_names,2,mkChar("Residuals"));
     SET_STRING_ELT(R_return_value_names,3,mkChar("StdErrors"));
     SET_STRING_ELT(R_return_value_names,4,mkChar("Scale"));
     setAttrib(R_return_value_cur, R_NamesSymbol,R_return_value_names);
     UNPROTECT(2);
   
     pthread_mutex_lock (&mutex_R);
     SET_VECTOR_ELT(*(args->R_return_value),j,R_return_value_cur);
     pthread_mutex_unlock (&mutex_R);
  }
}
#endif





SEXP R_sub_rcModelSummarize_plm(SEXP RMatrix, SEXP R_rowIndexList, SEXP PsiCode, SEXP PsiK, SEXP Scales){

  SEXP R_return_value;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer, *buffer2;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;
#ifdef USE_PTHREADS
  int t, returnCode, chunk_size, num_threads = 1;
  double chunk_size_d, chunk_tot_d;
  char *nthreads;
  pthread_attr_t attr;
  pthread_t *threads;
  struct loop_data *args;
  void *status; 
#ifdef PTHREAD_STACK_MIN
  size_t stacksize = PTHREAD_STACK_MIN + 0x4000;
#else
  size_t stacksize = 0x8000;
#endif
#else

  SEXP R_return_value_cur;

  SEXP R_weights;
  SEXP R_residuals;
  SEXP R_beta;
  SEXP R_SE;
  SEXP R_scale;

  SEXP R_return_value_names;

  double *Ymat;

  double *beta;
  double *residuals;
  double *weights;
  double *se;

  double scale=-1.0;
  double *scaleptr;

  double residSE;

#endif

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_return_value = allocVector(VECSXP,length_rowIndexList));
  
#ifdef  USE_PTHREADS
  nthreads = getenv(THREADS_ENV_VAR);
  if(nthreads != NULL){
    num_threads = atoi(nthreads);
    if(num_threads <= 0){
      error("The number of threads (enviroment variable %s) must be a positive integer, but the specified value was %s", THREADS_ENV_VAR, nthreads);
    }
  }
  threads = (pthread_t *) Calloc(num_threads, pthread_t);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setstacksize (&attr, stacksize);
  
  /* this code works out how many threads to use and allocates ranges of subColumns to each thread */
  /* The aim is to try to be as fair as possible in dividing up the matrix */
  /* A special cases to be aware of: 
     1) Number of subColumns is less than the number of threads
  */
  
  if (num_threads < length_rowIndexList){
    chunk_size = length_rowIndexList/num_threads;
    chunk_size_d = ((double) length_rowIndexList)/((double) num_threads);
  } else {
    chunk_size = 1;
    chunk_size_d = 1;
  }

  if(chunk_size == 0){
    chunk_size = 1;
  }
  args = (struct loop_data *) Calloc((length_rowIndexList < num_threads ? length_rowIndexList : num_threads), struct loop_data);

  args[0].matrix = matrix;
  args[0].R_return_value = &R_return_value;
  args[0].R_rowIndexList = &R_rowIndexList;
  args[0].PsiCode = &PsiCode;
  args[0].PsiK = &PsiK;
  args[0].Scales = &Scales;
  args[0].rows = rows;  
  args[0].cols = cols;
  args[0].length_rowIndexList = length_rowIndexList;

  pthread_mutex_init(&mutex_R, NULL);

  t = 0; /* t = number of actual threads doing work */
  chunk_tot_d = 0;
  for (i=0; floor(chunk_tot_d+0.00001) < length_rowIndexList; i+=chunk_size){
     if(t != 0){
       memcpy(&(args[t]), &(args[0]), sizeof(struct loop_data));
     }

     args[t].start_row = i;     
     /* take care of distribution of the remainder (when #chips%#threads != 0) */
     chunk_tot_d += chunk_size_d;
     // Add 0.00001 in case there was a rounding issue with the division
     if(i+chunk_size < floor(chunk_tot_d+0.00001)){
       args[t].end_row = i+chunk_size;
       i++;
     }
     else{
       args[t].end_row = i+chunk_size-1;
     }
     t++;
  }

  
  for (i =0; i < t; i++){
     returnCode = pthread_create(&threads[i], &attr, sub_rcModelSummarize_plm_group, (void *) &(args[i]));
     if (returnCode){
         error("ERROR; return code from pthread_create() is %d\n", returnCode);
     }
  }
  /* Wait for the other threads */
  for(i = 0; i < t; i++){
      returnCode = pthread_join(threads[i], &status);
      if (returnCode){
         error("ERROR; return code from pthread_join(thread #%d) is %d, exit status for thread was %d\n", 
               i, returnCode, *((int *) status));
      }
  }

  pthread_attr_destroy(&attr);  
  pthread_mutex_destroy(&mutex_R);
  Free(threads);
  Free(args);  
#else     

  for (j =0; j < length_rowIndexList; j++){    

    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));

    PROTECT(R_return_value_cur = allocVector(VECSXP,5));
    PROTECT(R_beta = allocVector(REALSXP, ncur_rows + cols));
    PROTECT(R_weights = allocMatrix(REALSXP,ncur_rows,cols));
    PROTECT(R_residuals = allocMatrix(REALSXP,ncur_rows,cols));
    PROTECT(R_SE = allocVector(REALSXP,ncur_rows+cols)); 
    PROTECT(R_scale = allocVector(REALSXP,1));

    SET_VECTOR_ELT(R_return_value_cur,0,R_beta);
    SET_VECTOR_ELT(R_return_value_cur,1,R_weights);
    SET_VECTOR_ELT(R_return_value_cur,2,R_residuals);
    SET_VECTOR_ELT(R_return_value_cur,3,R_SE);
    SET_VECTOR_ELT(R_return_value_cur,4,R_scale);

    UNPROTECT(5);

    beta = NUMERIC_POINTER(R_beta);
    residuals = NUMERIC_POINTER(R_residuals);
    weights = NUMERIC_POINTER(R_weights);
    se = NUMERIC_POINTER(R_SE);
    

   scaleptr = NUMERIC_POINTER(R_scale);


    if (isNull(Scales)){
      scaleptr[0] = -1.0;
    } else if (length(Scales) != cols) {
      scaleptr[0] = NUMERIC_POINTER(Scales)[0];
    }


    Ymat = Calloc(ncur_rows*cols,double);
    
    
    for (k = 0; k < cols; k++){
        for (i =0; i < ncur_rows; i++){
     	    Ymat[k*ncur_rows + i] = args->matrix[k*args->rows + cur_rows[i]];  
        }
    } 

    rlm_fit_anova_scale(Ymat, ncur_rows, cols, scaleptr, beta, residuals, weights, PsiFunc(asInteger(PsiCode)),asReal(PsiK), 20, 0);
  
    rlm_compute_se_anova(Ymat, ncur_rows, cols, beta, residuals, weights,se, (double *)NULL, &residSE, 4, PsiFunc(asInteger(PsiCode)),asReal(PsiK));


  

     beta[ncur_rows+cols -1] = 0.0;
  

     for (i = cols; i < ncur_rows + cols -1; i++)
        beta[ncur_rows+cols -1]-=beta[i];

     Free(Ymat);
     PROTECT(R_return_value_names= allocVector(STRSXP,5));
     SET_STRING_ELT(R_return_value_names,0,mkChar("Estimates"));
     SET_STRING_ELT(R_return_value_names,1,mkChar("Weights"));
     SET_STRING_ELT(R_return_value_names,2,mkChar("Residuals"));
     SET_STRING_ELT(R_return_value_names,3,mkChar("StdErrors"));
     SET_STRING_ELT(R_return_value_names,4,mkChar("Scale"));
     setAttrib(R_return_value_cur, R_NamesSymbol,R_return_value_names);
     UNPROTECT(2);
     SET_VECTOR_ELT(R_return_value,j,R_return_value_cur);
  }
#endif
  UNPROTECT(1);
  return R_return_value;
}






