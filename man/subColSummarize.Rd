\name{subColSummarize}
\alias{subColSummarizeAvg}
\alias{subColSummarizeAvgLog}
\alias{subColSummarizeBiweight}
\alias{subColSummarizeBiweightLog}
\alias{subColSummarizeLogAvg}
\alias{subColSummarizeLogMedian}
\alias{subColSummarizeMedian}
\alias{subColSummarizeMedianLog}
\alias{subColSummarizeMedianpolish}
\alias{subColSummarizeMedianpolishLog}
\alias{convert.group.labels}
\title{Summarize matrix columns by row groups}
\description{These functions partition the rows of a matrix into groups 
  (specified via \code{group.labels}) and then summarize each column independently 
  within each group. They use the same summarization methods as \code{\link{colSummarize}} 
  (mean, median, Tukey Biweight, median polish, and log-transformed variants) but apply 
  them to subsets of rows. The result is a matrix where each row represents the summarized 
  values for one group across all columns.
}
\usage{subColSummarizeAvg(y, group.labels)
       subColSummarizeAvgLog(y, group.labels)
       subColSummarizeBiweight(y, group.labels)
       subColSummarizeBiweightLog(y, group.labels)
       subColSummarizeLogAvg(y, group.labels)
       subColSummarizeLogMedian(y, group.labels)
       subColSummarizeMedian(y, group.labels)
       subColSummarizeMedianLog(y, group.labels)
       subColSummarizeMedianpolish(y, group.labels)
       subColSummarizeMedianpolishLog(y, group.labels)
       convert.group.labels(group.labels)
}
\arguments{
  \item{y}{A numeric matrix where each column represents an independent group or sample 
    and each row represents an observation or measurement within that group. 
    The matrix should contain at least 1 row and 1 column.}
  \item{group.labels}{A vector of length \code{nrow(y)} that assigns each row to a group. 
    Can be any atomic vector (numeric, character, logical) that will be converted to a factor.
    Rows with \code{NA} values in \code{group.labels} will be excluded from summarization.
    The names of the resulting groups will be taken from the factor level names.}
}
\value{
  A numeric matrix with dimensions: number of groups × number of columns.
  Each row corresponds to one group (with row names from the factor levels of \code{group.labels}),
  and each column corresponds to a column from the input matrix.
  The matrix contains the summarized values for each group and column combination.
}
\details{
  These functions are designed to summarize the columns of a matrix
  where the rows of the matrix are assigned to groups. The summarization
  is by column across all rows in each group.
  \itemize{
    \item{subColSummarizeAvg}{Summarize by taking mean}
    \item{subColSummarizeAvgLog}{\code{log2} transform the data and
      then take means in column-wise manner}
    \item{subColSummarizeBiweight}{Use a one-step Tukey Biweight to
      summarize columns}
    \item{subColSummarizeBiweightLog}{\code{log2} transform the data and
      then use a one-step Tukey Biweight to
      summarize columns}
    \item{subColSummarizeLogAvg}{Summarize by taking mean and then
      taking \code{log2}}
    \item{subColSummarizeLogMedian}{Summarize by taking median and then
      taking \code{log2}}
    \item{subColSummarizeMedian}{Summarize by taking median}
    \item{subColSummarizeMedianLog}{\code{log2} transform the data and
      then summarize by taking median}
    \item{subColSummarizeMedianpolish}{Use the median polish to
      summarize each column, by also using a row effect (not returned)}
    \item{subColSummarizeMedianpolishLog}{\code{log2} transform the
      data and then use the median polish to summarize each column, by
      also using a row effect (not returned)}
  }
  
}
\section{When to Use subColSummarize}{

  Use \code{subColSummarize*} functions when your data has a natural grouping 
  structure within the rows, and you want to compute summary statistics separately 
  for each group. Common use cases include:
  
  \describe{
    \item{Stratified analysis}{When you have multiple samples/replicates within each 
      treatment group or experimental condition, use these functions to summarize 
      by condition rather than across all rows at once.}
    
    \item{Quality control grouping}{When rows represent different probe sets, 
      detection methods, or measurement batches, group them to compute batch-specific 
      summaries and identify systematic differences between groups.}
    
    \item{Multi-level experiments}{When your data represents a hierarchical structure 
      (e.g., multiple technical replicates per biological sample), group by the 
      higher level (biological sample) to obtain summary values accounting for 
      within-group variation.}
    
    \item{Robust comparison}{Use with robust methods (Biweight, median polish) 
      to compare groups while being resistant to outliers within each group.}
  }
  
  \strong{Example:} If you have expression data with 3 probes × 4 samples where 
  the probes are grouped as \{probe1, probe2, probe1, probe2\}, you could group 
  by probe to get summary expression values per probe per sample, effectively 
  combining replicates from the same probe.
}
\examples{
# Set random seed for reproducible examples
set.seed(1234)

# Create sample data: 20 rows (10 in group 1, 10 in group 2) and 5 columns
y <- matrix(c(10 + rnorm(50), 20 + rnorm(50)), 20, 5, byrow = TRUE)
groups <- c(rep(1, 10), rep(2, 10))

# Simple methods - returns a matrix with one row per group
subColSummarizeAvg(y, groups)
subColSummarizeMedian(y, groups)

# Log-transformed average methods
subColSummarizeAvgLog(y, groups)
subColSummarizeLogAvg(y, groups)

# Log-transformed median methods
subColSummarizeLogMedian(y, groups)
subColSummarizeMedianLog(y, groups)

# Robust methods - resistant to outliers within groups
subColSummarizeBiweight(y, groups)
subColSummarizeBiweightLog(y, groups)

# Median polish methods - accounts for group and column effects
subColSummarizeMedianpolish(y, groups)
subColSummarizeMedianpolishLog(y, groups)

# Examine the output: a matrix with 2 rows (one per group) and 5 columns
result <- subColSummarizeAvg(y, groups)
dim(result)           # 2 groups × 5 columns
rownames(result)      # Group names
}
\seealso{
  \code{\link{colSummarize}} for summarizing entire columns without grouping
}
\author{B. M. Bolstad <bmb@bmbolstad.com>}
\keyword{univar}