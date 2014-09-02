/************************************************************************
 **
 ** median.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Feb 6, 2003  (but based on earlier work from median_log.c)
 **
 ** Copyright (C) 2007   Ben Bolstad
 **
 ** last modified: Sep 16, 2007
 **
 ** License: LGPL V2 (same as the rest of the preprocessCore package)
 **
 ** General discussion
 **
 ** Implement median log2 pm summarization.
 **
 ** Sep 16, 2007 - initial version
 ** Sep, 2014 - Change to size_t where appropriate. Code documentation cleanup
 **
 ************************************************************************/

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

#include "rma_common.h"
#include "median.h"


/***************************************************************************
 **
 ** double colmedian_wrapper(double *x, size_t length)
 **
 ** double *x - a vector of PM intensities 
 ** size_t length - length of *x
 **
 ** take the median of PM intensities.
 **
 ***************************************************************************/

static double colmedian_wrapper(double *x, size_t length){

  double med = 0.0;
  
  med = median_nocopy(x,length);
 
  return (med);    
}

/***************************************************************************
 **
 ** double ColMedian(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute Median expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 ** double *resultsSE - already allocated location to store expression measures SE estimates (cols length)
 **
 ***************************************************************************/

/*! \brief Compute the median and SE of the median for subset of rows
 * 
 *  Given a data matrix of probe intensities compute median and SE of this estimate
 *  on a column by column basis using only a specified subset of rows. 
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows indices specifying which rows in the matrix to use
 * @param results pre-allocated space to store output medians. Should be of length cols
 * @param nprobes the number of elements in cur_rows
 * @param resultsSE pre-allocated space to store SE of medians. Should be of length cols
 *
 *  
 */

void ColMedian(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE){

  size_t i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = colmedian_wrapper(&z[j*nprobes],nprobes); 
    resultsSE[j] = R_NaReal;
  }
  Free(z);
}



/***************************************************************************
 **
 ** double ColMedian_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute Median expression measure. 
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

/*! \brief Compute the median 
 * 
 *  Given a data matrix of probe intensities compute median measure 
 *  on a column by column basis. 
 *  
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

void ColMedian_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes){

  size_t i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = colmedian_wrapper(&z[j*nprobes],nprobes);

  }
  Free(z);
}



/*! \brief Compute the median and SE of the median
 * 
 *  Given a data matrix of probe intensities compute median measure and SE of this estimate
 *  on a column by column basis. The sample standard error is also computed. l
 *  
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

void colmedian(double *data, size_t rows, size_t cols, double *results, double *resultsSE){

  size_t i,j;
  double *buffer = Calloc(rows, double);
  
  for (j=0; j < cols; j++){
    for (i = 0; i < rows; i++){
      buffer[i] = data[j*rows + i];
    }
    results[j] = colmedian_wrapper(buffer,rows); 
    resultsSE[j] = R_NaReal;
  }

  Free(buffer);
}



/*! \brief Compute the median and SE of the median
 * 
 *  Given a data matrix of probe intensities compute median measure and SE of this estimate
 *  on a column by column basis. The sample standard error is also computed. On output the data matrix will
 *  be changed.
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

void colmedian_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE){

  size_t i,j;
    
  for (j=0; j < cols; j++){
    results[j] = colmedian_wrapper(&data[j*rows],rows); 
    resultsSE[j] = R_NaReal;
  }

}
