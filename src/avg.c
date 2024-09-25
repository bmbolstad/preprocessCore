/************************************************************************
 **
 ** avg.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Sep 16, 2007  (but based on earlier work from Nov avg_log.c)
 **
 ** Copyright (C) 2007-2014 Ben Bolstad
 **
 ** last modified: Sept 1, 2014
 **
 ** License: LGPL V2 (same as the rest of the preprocessCore package)
 **
 ** General discussion
 **
 ** Implement average summarization
 **
 ** This file provides functions that take the mean (arithmetic average) by
 ** column of matrixes. They may also return the standard error estimate of the mean.
 ** 
 ** There are four main functions (that are exposed to outside this file):
 ** colaverage  - computes averages of each column
 ** colaverage_no_copy - computes averages of each column (does not allocate extra space, which means may change values in input matrix)
 ** ColAverage  - computes averages (and SE of average) of each column using only a subset of rows (with subset specified and identical across columns)
 ** ColAverage_noSE - computes averages of each column using only a subset of rows (with subset specified and identical across columns)
 **
 ** History
 ** Sep 16, 2007 - Initial version
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

#include "avg.h"



/***************************************************************************
 **
 ** double Avg(double *x, size_t length)
 **
 ** double *x - a vector of PM intensities  
 ** size_t length - length of *x
 **
 ** take the average of input intensities.
 **
 ***************************************************************************/

static double Avg(double *x, size_t length){
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
 ** static double AvgSE(double *x, double mean, size_t length)
 **
 ** double *x - a vector of PM intensities 
 ** double mean - the mean of x computed using Avg above
 ** int length - length of *x
 **
 ** compute the standard error of the average of log2 PM intensities.
 ** 
 **
 ***************************************************************************/

static double AvgSE(double *x, double mean, size_t length){
  int i;
  double sum = 0.0;

  for (i=0; i < length; i++){
    sum = sum + (x[i]- mean)*(x[i] - mean);
  }
  
  sum = sqrt(sum/(double)(length -1));
  sum = sum/sqrt((double)length);

  return (sum);    
}



/***************************************************************************
 ** 
 ** void colaverage_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, compute averages in column wise manner. also return SE of mean 
 **      
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

/*! \brief Compute the mean and SE of the mean
 * 
 *  Given a data matrix of probe intensities compute average expression measure and SE of this estimate
 *  on a column by column basis. Specifically, the arithmetic mean
 *  is computed for each column. The sample standard error is also computed. This function guarantees that 
 *  no additional memory is temporarily allocated to copy the input data matrix. However, this means that
 *  on output the input matrix may be changed.
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

void colaverage_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE){
  int j;

  for (j = 0; j < cols; j++){
    results[j] = Avg(&data[j*rows],rows);
    resultsSE[j] = AvgSE(&data[j*rows],results[j],rows);
  } 
}



/***************************************************************************
 ** 
 ** void colaverage(double *data, size_t rows, size_t cols, double *results, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, compute averages in column wise manner 
 **      
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

/*! \brief Compute the mean and SE of the mean
 * 
 *  Given a data matrix of probe intensities compute average expression measure and SE of this estimate
 *  on a column by column basis. Specifically, the arithmetic mean
 *  is computed for each column. The sample standard error is also computed. On output the data matrix will
 *  be unchanged.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param results pre-allocated space to store output averages. Should be of length cols
 * @param resultsSE pre-allocated space to store SE of averages. Should be of length cols
 *
 *  
 */

void colaverage(double *data, size_t rows, size_t cols, double *results, double *resultsSE){
  int i,j;
  double *z = R_Calloc(rows,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < rows; i++){
      z[i] = data[j*rows + i];  
    }
    results[j] = Avg(z,rows);
    resultsSE[j] = AvgSE(z,results[j],rows);
  } 
  R_Free(z);

}



/***************************************************************************
 **
 ** double ColAverage(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset. 
 ** double *resultsSE - already allocated location to store expression measures SE (cols length)
 **
 ***************************************************************************/

/*! \brief Given a data matrix of probe intensities, and a list of rows in the matrix 
 *      corresponding to a single probeset, compute average expression measure. 
 *      Note that data is a probes by chips matrix. Also compute SE estimates
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows a vector containing row indices to use
 * @param results pre-allocated space to store output averages. Should be of length cols
 * @param nprobes number of probes in current set
 * @param resultsSE pre-allocated space to store SE of averages. Should be of length cols
 *
 *  
 */

void ColAverage(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE){
  int i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = Avg(&z[j*nprobes],nprobes);
    resultsSE[j] = AvgSE(&z[j*nprobes],results[j],nprobes);
  }

  R_Free(z);
}



/***************************************************************************
 **
 ** void ColAverage_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
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
 *      corresponding to a single probeset, compute average expression measure. 
 *      Note that data is a probes by chips matrix. 
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows a vector containing row indices to use
 * @param results pre-allocated space to store output averages. Should be of length cols
 * @param nprobes number of probes in current set
 *
 *  
 */

void ColAverage_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes){
  int i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = Avg(&z[j*nprobes],nprobes);
  }
  R_Free(z);
}
