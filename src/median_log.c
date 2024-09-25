/************************************************************************
 **
 ** median_logPM.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Feb 6, 2003  (but based on earlier work from Nov 2002)
 **
 ** Copyright (C) 2003-2014   Ben Bolstad
 **
 ** last modified: Feb 6, 2003
 **
 ** License: LGPL V2 (same as the rest of the preprocessCore package)
 **
 ** General discussion
 **
 ** Implement median log2 pm summarization.
 **
 ** Feb 6, 2003 - Initial version of this summarization method
 ** Feb 24, 2003 - Remove unused variable in i from MedianLog
 ** Jul 23, 2003 - add SE parameter (but not yet implemented)
 ** Oct 10, 2003 - added PLM version
 ** Sept 9, 2007 - branch out of affyPLM into a new package preprocessCore
 ** Sept, 2014 - Documentation clean up. Change to size_t where appropriate
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

#include "rma_common.h"
#include "median_log.h"

/***************************************************************************
 **
 ** double MedianLog(double *x, size_t length)
 **
 ** double *x - a vector of PM intensities 
 ** int length - length of *x
 **
 ** take the log2 of the median of PM intensities.
 **
 ***************************************************************************/

static double median_log(double *x, size_t length){

  double med = 0.0;
  
  med = median_nocopy(x,length);
 
  return (med);    
}



/***************************************************************************
 **
 ** double MedianLogPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute log2 Median expression measure. 
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

/*! \brief  \f$\log_2\f$ transform the data and compute the median 
 * 
 *  Given a data matrix of probe intensities \f$\log_2\f$ transform it and then compute the median. Also compute SE of this estimate
 *  on a column by column basis using only a specified subset of rows. Specifically, the median of each column is based on
 *  \f$\log_2\f$ transformed data. The sample standard error is also computed. 
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows indices specifying which rows in the matrix to use
 * @param results pre-allocated space to store output log2 medians. Should be of length cols
 * @param nprobes the number of rows in cur_rows
 * @param resultsSE pre-allocated space to store SE of log2 medians. Should be of length cols
 *
 *  
 */

void MedianLog(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes, double *resultsSE){

  size_t i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = median_log(&z[j*nprobes],nprobes); 
    resultsSE[j] = R_NaReal;
  }
  R_Free(z);
}



/***************************************************************************
 **
 ** double MedianLogPM_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute log2 Median expression measure. 
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

/*! \brief  \f$\log_2\f$ transform the data and compute the median 
 * 
 *  Given a data matrix of probe intensities \f$\log_2\f$ transform it and then compute the median on a column by column basis using only a specified subset of rows. 
 *  Specifically, the median of each column is based on \f$\log_2\f$ transformed data.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param cur_rows indices specifying which rows in the matrix to use
 * @param results pre-allocated space to store output log2 medians. Should be of length cols
 * @param nprobes the number of rows in cur_rows
 *  
 */

void MedianLog_noSE(double *data, size_t rows, size_t cols, int *cur_rows, double *results, size_t nprobes){

  size_t i,j;
  double *z = R_Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = median_log(&z[j*nprobes],nprobes);

  }
  R_Free(z);
}



/*! \brief compute the median for each column of \f$\log_2\f$ transformed data.
 * 
 *  Given a data matrix of probe intensities \f$\log_2\f$ transform it then compute median of each column. Also produce the SE of this estimate
 *  on a column by column basis. Specifically, the median is computed for each column and then \f$\log_2\f$ transformed.
 *  The sample standard error is also computed. On output the data matrix will
 *  be unchanged.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param results pre-allocated space to store output log2 medians. Should be of length cols
 * @param resultsSE pre-allocated space to store SE of log2 medians. Should be of length cols
 *
 *  
 */

void medianlog(double *data, size_t rows, size_t cols, double *results, double *resultsSE){

  size_t i,j;
  double *buffer = R_Calloc(rows, double);
  
  for (j=0; j < cols; j++){
    for (i = 0; i < rows; i++){
      buffer[i] = log(data[j*rows + i])/log(2.0);
    }
    results[j] = median_log(buffer,rows); 
    resultsSE[j] = R_NaReal;
  }

  R_Free(buffer);
}



/*! \brief compute the median for each column of \f$\log_2\f$ transformed data.
 * 
 *  Given a data matrix of probe intensities \f$\log_2\f$ transform it then compute median of each column. Also produce the SE of this estimate
 *  on a column by column basis. Specifically, the median is computed for each column and then \f$\log_2\f$ transformed.
 *  The sample standard error is also computed. On output the data matrix will
 *  be changed.
 *    
 *
 * @param data a matrix containing data stored column-wise stored in rows*cols length of memory
 * @param rows the number of rows in the matrix 
 * @param cols the number of columns in the matrix
 * @param results pre-allocated space to store output log2 medians. Should be of length cols
 * @param resultsSE pre-allocated space to store SE of log2 medians. Should be of length cols
 *
 *  
 */

void medianlog_no_copy(double *data, size_t rows, size_t cols, double *results, double *resultsSE){

  size_t i,j;
    
  for (j=0; j < cols; j++){
    for (i = 0; i < rows; i++){
      data[j*rows + i]= log(data[j*rows + i])/log(2.0);
    }
    results[j] = median_log(&data[j*rows],rows); 
    resultsSE[j] = R_NaReal;
  }

}
