\name{colSummarize}
\alias{colSummarizeAvg}
\alias{colSummarizeAvgLog}
\alias{colSummarizeBiweight}
\alias{colSummarizeBiweightLog}
\alias{colSummarizeLogAvg}
\alias{colSummarizeLogMedian}
\alias{colSummarizeMedian}
\alias{colSummarizeMedianLog}
\alias{colSummarizeMedianpolish}
\alias{colSummarizeMedianpolishLog}

\title{Summarize the columns of matrices}
\description{Compute column-wise summary values of a matrix.
}
\usage{
colSummarizeAvg(y)
colSummarizeAvgLog(y)
colSummarizeBiweight(y)
colSummarizeBiweightLog(y)
colSummarizeLogAvg(y)
colSummarizeLogMedian(y)
colSummarizeMedian(y)
colSummarizeMedianLog(y)
colSummarizeMedianpolish(y)
colSummarizeMedianpolishLog(y)

}
\arguments{
  \item{y}{A numeric matrix. Each column represents an independent group or sample 
    and each row represents an observation or measurement within that group. 
    The matrix should contain at least 1 row and 1 column. Values must be numeric; 
    \code{NA} and \code{NaN} values will be propagated to the results.}
}
\value{
  A named list with two elements:
  \item{Estimates}{A numeric vector of length \code{ncol(y)} containing 
    the summary value for each column.}
  \item{StdErrors}{A numeric vector of length \code{ncol(y)} containing 
    the standard error estimate for each column. Note that 
    \code{colSummarizeLogMedian} and \code{colSummarizeMedianLog} return 
    \code{NA} for all elements of this vector.}
}
\details{This group of functions summarize the columns of a given
  matrix. Note that \code{colSummarizeLogMedian} and \code{colSummarizeMedianLog} 
  return \code{NA} for standard errors, as SE is not computed for these methods.

  \itemize{
    \item{\code{colSummarizeAvg}}{Take means in column-wise manner}
    \item{\code{colSummarizeAvgLog}}{\code{log2} transform the data and
      then take means in column-wise manner}
    \item{\code{colSummarizeBiweight}}{Summarize each column using a one
      step Tukey Biweight procedure}
    \item{\code{colSummarizeBiweightLog}}{\code{log2} transform the data
      and then summarize each column using a one step Tukey Biweight
      procedure}
    \item{\code{colSummarizeLogAvg}}{Compute the mean of each column and
      then \code{log2} transform it}
    \item{\code{colSummarizeLogMedian}}{Compute the median of each
      column and then \code{log2} transform it}
    \item{\code{colSummarizeMedian}}{Compute the median of each column}
    \item{\code{colSummarizeMedianLog}}{\code{log2} transform the data
      and then summarize each column using the median}
    \item{\code{colSummarizeMedianpolish}}{Use the median polish to
      summarize each column, by also using a row effect (not returned)}
    \item{\code{colSummarizeMedianpolishLog}}{\code{log2} transform the
      data and then use the median polish to summarize each column, by
      also using a row effect (not returned)}
  }
}
\section{Choosing a Summarization Method}{

  The choice of method depends on your data characteristics and analysis goals:

  \describe{
    \item{Simple methods (Avg, Median)}{
      Use for clean data without significant outliers. 
      \code{colSummarizeAvg} is suitable for approximately normally distributed data, 
      while \code{colSummarizeMedian} is more robust to skewed distributions and outliers.
      Both return results on the original data scale.
    }
    
    \item{Log-transformed methods (AvgLog, LogAvg, MedianLog, LogMedian)}{
      Use when data exhibits multiplicative noise or exponential behavior (common in 
      biological measurements). \code{colSummarizeAvgLog} and \code{colSummarizeLogAvg} 
      are mathematically related through logarithm properties; \code{colSummarizeAvgLog} 
      is often preferred for right-skewed data. Note that \code{colSummarizeLogMedian} 
      and \code{colSummarizeMedianLog} do not return standard error estimates.
    }
    
    \item{Robust methods (Biweight, Biweight with log)}{
      Use when you expect outliers or want resistance to extreme values. 
      The Tukey Biweight method provides robust estimation with less influence from outliers 
      compared to means. Choose \code{colSummarizeBiweightLog} if your data is best 
      modeled on a log scale.
    }
    
    \item{Median polish methods (Medianpolish, Medianpolish with log)}{
      Use when you have a two-way structure (both row and column effects) and want 
      to account for systematic differences between rows. This method simultaneously 
      estimates column effects (returned as Estimates) while accounting for row effects.
    }
  }
}
\examples{
# Set random seed for reproducible examples
set.seed(1234)

# Create sample data: 20 observations per column, 5 columns
y <- matrix(10 + rnorm(100), 20, 5)

# Simple methods - returns mean and standard error for each column
colSummarizeAvg(y)
colSummarizeMedian(y)

# Log-transformed average methods (useful for right-skewed data)
colSummarizeAvgLog(y)
colSummarizeLogAvg(y)

# Log-transformed median methods (no standard errors returned)
colSummarizeLogMedian(y)
colSummarizeMedianLog(y)

# Robust methods - resistant to outliers using Tukey Biweight
colSummarizeBiweight(y)
colSummarizeBiweightLog(y)

# Median polish methods - accounts for both column and row effects
colSummarizeMedianpolish(y)
colSummarizeMedianpolishLog(y)

# Example of extracting results
result <- colSummarizeAvg(y)
head(result$Estimates)    # Column summaries
head(result$StdErrors)    # Standard errors
}
\seealso{
  \code{\link{subColSummarize}} for summarizing subsets of columns
}
\author{B. M. Bolstad   \email{bmb@bmbolstad.com}}
\keyword{univar}