/*********************************************************************
 **
 ** file: R_subColSummarize.c
 **
 ** Aim: Code which provides .Call() interfaces to the subcolumn 
 ** summarization code.
 **
 ** Copyright (C) 2007 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bmb@bmbolstad.com>
 ** 
 ** created on: Sep 15, 2007
 **
 ** History
 ** Sep 18, 2007 - Initial version
 **
 **
 *********************************************************************/

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>


#include "avg_log.h"
#include "log_avg.h"
#include "avg.h"

#include "biweight.h"

#include "median_log.h"
#include "log_median.h"
#include "median.h"

#include "medianpolish.h"



SEXP R_subColSummarize_avg_log(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    AverageLog_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}


SEXP R_subColSummarize_log_avg(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    LogAverage_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}


SEXP R_subColSummarize_avg(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    ColAverage_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}






SEXP R_subColSummarize_biweight_log(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    TukeyBiweight_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}




SEXP R_subColSummarize_biweight(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    TukeyBiweight_no_log_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}






SEXP R_subColSummarize_median_log(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    MedianLog_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}





SEXP R_subColSummarize_log_median(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    LogMedian_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}




SEXP R_subColSummarize_median(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    ColMedian_noSE(matrix, rows, cols, cur_rows, buffer, ncur_rows);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
  
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}




SEXP R_subColSummarize_medianpolish_log(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer, *buffer2;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);
  buffer2 = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    MedianPolish(matrix, rows, cols, cur_rows, buffer, ncur_rows, buffer2);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
   Free(buffer2);
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}





SEXP R_subColSummarize_medianpolish(SEXP RMatrix, SEXP R_rowIndexList){

  SEXP R_summaries;  
  SEXP dim1;

  double *matrix=NUMERIC_POINTER(RMatrix);
  double *results, *buffer, *buffer2;
  
  int *cur_rows;

  int rows, cols;
  int length_rowIndexList = LENGTH(R_rowIndexList);
  int ncur_rows;

  int i,j;

  PROTECT(dim1 = getAttrib(RMatrix,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  PROTECT(R_summaries = allocMatrix(REALSXP,length_rowIndexList,cols));
 
  results = NUMERIC_POINTER(R_summaries);
 
  buffer = Calloc(cols,double);
  buffer2 = Calloc(cols,double);

  for (j =0; j < length_rowIndexList; j++){    
    ncur_rows = LENGTH(VECTOR_ELT(R_rowIndexList,j)); 
    cur_rows = INTEGER_POINTER(VECTOR_ELT(R_rowIndexList,j));
    MedianPolish_no_log(matrix, rows, cols, cur_rows, buffer, ncur_rows, buffer2);
    
    for (i = 0; i < cols; i++){
      results[i*length_rowIndexList + j] = buffer[i];
    }
  }
   Free(buffer2);
  Free(buffer);
  UNPROTECT(1);
  return R_summaries;
}
