/************************************************************************
 **
 ** avg_log.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Jan 7, 2003  (but based on earlier work from Nov 2002)
 **
 ** Copyright (C) 2002-2014 Ben Bolstad
 **
 ** last modified: Jan 7, 2003
 **
 ** License: LGPL V2 (same as the rest of the preprocessCore package)
 **
 ** General discussion
 **
 ** Implement average log2 pm summarization, with or without normalization
 **
 **  
 ** There are four main functions (that are exposed to outside this file):
 ** averagelog  - computes averages of log2 of each column
 ** averagelog_no_copy - computes averages of log2 of each column (does not allocate extra space, which means may change values in input matrix)
 ** AverageLog  - computes averages (and SE of average) of log2 of each column using only a subset of rows (with subset specified and identical across columns)
 ** AverageLog_noSE - computes averages of log2 of each column using only a subset of rows (with subset specified and identical across columns)
 **
 ** Nov 2, 2002 - modify so that it will work efficently with affy2
 ** Jan 3, 2003 - Clean up commenting, prepare for integration in AffyExtensions
 ** Jan 7, 2003 - Make function standalone, to prepare for later combination into
 **               a more general framework.
 ** Jul 23, 2003 - add parameter for computing SE and SE implemented
 ** Oct 5, 2003 - add output_param
 ** Oct 10, 2003 - added threestepPLM version of this summary.
 ** May 19, 2007 - branch out of affyPLM into a new package preprocessCore, then restructure the code. Add doxygen style documentation
 ** May 26, 2007 - fix memory leak in average_log. add additional interfaces
 ** Sep 16, 2007 - fix error in how StdError is computed
 ** Sep 2014 - Change to size_t rather than int for variables indexing pointers. Improve code documentation.
 **
 **
 ************************************************************************/

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#include "avg_log.h"


/***************************************************************************
 **
 ** double AvgLog(double *x, size_t length)
 **
 ** double *x - a vector of PM intensities  (previously log2 transformed)
 ** size_t length - length of *x
 **
 ** take the average of log2 PM intensities.
 **
 ***************************************************************************/

static double AvgLog(double *x, size_t length){
  int i;
  double sum = 0.0;
  double mean = 0.0;

  for (i=0; i < length; i++){
    sum = sum + x[i];
  }
  
  mean = sum/(double)length;

  return (mean);    
}



/***************************************************************************
 **
 ** static double AvgLogSE(double *x, size_t length)
 **
 ** double *x - a vector of PM intensities (previously log2 transformed)
 ** double mean - the mean of x computed using AvgLog above
 ** size_t length - length of *x
 **
 ** compute the standard error of the average of log2 PM intensities.
 ** 
 **
 ***************************************************************************/

static double AvgLogSE(double *x, double mean, size_t length){
  int i;
  double sum = 0.0;

  for (i=0; i < length; i++){
    sum = sum + (x[i]- mean)*(x[i] - mean);
  }
  
  sum = sqrt(sum/(double)(length-1));
  sum = sum/sqrt((double)length);

  return (sum);    
}



/***************************************************************************
 ** 
 ** void averagelog_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, compute average of log2 values in column wise manner 
 **      no additional memory allocation is done, input matrix may be changed on output
 **
 ** double *data - Probe intensity matrix
 ** size_t rows - number of rows in matrix *data (probes)
 ** size_t cols - number of cols in matrix *data (chips)
 ** double *results - already allocated location to store expression measures (cols length)
 ** double *resultsSE - already allocated location to store expression measures standard error (cols length)
 **
 ***************************************************************************/

/*! \brief log2 transform and then compute the mean and SE of the mean
 * 
 *  Given a data matrix of probe intensities compute average log2 expression measure and SE of this estimate
 *  on a column by column basis. Specifically, each element is log2 transformed, then the arithmetic mean
 *  is computed for each column. The sample standard error is also computed. This function guarantees that 
 *  no additional memory is temporarily allocated to copy the input data matrix. However, this means that
 *  on output the input matrix will be changed.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param results pre-allocated space to store output log2 averages. Should be of length cols
 * @param resultsSE pre-allocated space to store SE of log2 averages. Should be of length cols
 *
 *  
 */

void averagelog_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE){
  int i,j;

  for (j = 0; j < cols; j++){
    for (i =0; i < rows; i++){
      data[j*rows + i] = log(data[j*rows + i])/log(2.0);  
    }
    results[j] = AvgLog(&data[j*rows],rows);
    resultsSE[j] = AvgLogSE(&data[j*rows],results[j],rows);
  } 

}



/***************************************************************************
 ** 
 ** void averagelog(double *data, size_t rows, size_t cols, double *results, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, compute average of log2 values in column wise manner 
 **      
 **
 ** double *data - Probe intensity matrix
 ** size_t rows - number of rows in matrix *data (probes)
 ** size_t cols - number of cols in matrix *data (chips)
 ** double *results - already allocated location to store expression measures (cols length)
 ** double *resultsSE - already allocated location to store expression measures standard error (cols length)
 **
 ***************************************************************************/

/*! \brief log2 transform and then compute the mean and SE of the mean
 * 
 *  Given a data matrix of probe intensities compute average log2 expression measure and SE of this estimate
 *  on a column by column basis. Specifically, each element is log2 transformed, then the arithmetic mean
 *  is computed for each column. The sample standard error is also computed. On output the data matrix will
 *  be unchanged.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param results pre-allocated space to store output log2 averages. Should be of length cols
 * @param resultsSE pre-allocated space to store SE of log2 averages. Should be of length cols
 *
 *  
 */

void averagelog(double *data, size_t rows, size_t cols, double *results, double *resultsSE){
  int i,j;
  double *z = R_Calloc(rows,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < rows; i++){
      z[i] = log(data[j*rows + i])/log(2.0);  
    }
    results[j] = AvgLog(z,rows);
    resultsSE[j] = AvgLogSE(z,results[j],rows);
  } 
  R_Free(z);

}



/***************************************************************************
 **
 ** double AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average log2 expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 ** double *resultsSE - already allocated location to store expression measures standard errors (cols length)
 **
 ***************************************************************************/

/*! \brief Given a data matrix of probe intensities, and a list of rows in the matrix 
 *      corresponding to a single probeset, compute average log2 expression measure. 
 *      Note that data is a probes by chips matrix. Also compute SE estimates
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows a vector containing row indices to use
 * @param results pre-allocated space to store output log2 averages. Should be of length cols
 * @param nprobes number of probes in current set
 * @param resultsSE pre-allocated space to store SE of log2 averages. Should be of length cols
 *
 *  
 */

void AverageLog(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE){
  int i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = AvgLog(&z[j*nprobes],nprobes);
    resultsSE[j] = AvgLogSE(&z[j*nprobes],results[j],nprobes);
  }

  R_Free(z);
}



/***************************************************************************
 **
 ** double AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average log2 expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

/*! \brief Given a data matrix of probe intensities, and a list of rows in the matrix 
 *      corresponding to a single probeset, compute average log2 expression measure. 
 *      Note that data is a probes by chips matrix. 
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows a vector containing row indices to use
 * @param results pre-allocated space to store output log2 averages. Should be of length cols
 * @param nprobes number of probes in current set
 *
 *  
 */

void AverageLog_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes){
  int i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = AvgLog(&z[j*nprobes],nprobes);
  }
  R_Free(z);
}
