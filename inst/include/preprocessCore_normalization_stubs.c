#include <Rinternals.h>
#include <R_ext/Rdynload.h>



/*! \brief Quantile normalize the columns of a matrix
 *
 *  
 * @param data a matrix to be quantile normalized. On exit will be normalized
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 *
 */

int qnorm_c(double *data, int *rows, int *cols){


  static int(*fun)(double*, int*, int*) = NULL;
  
  if (fun == NULL)
    fun =  (int(*)(double*, int*, int*))R_GetCCallable("preprocessCore","qnorm_c");
  
  return fun(data,rows,cols);

}
