\name{normalize.quantiles.robust}
\alias{normalize.quantiles.robust}
\alias{normalize.AffyBatch.quantiles.robust}
\title{Robust Quantile Normalization}
\description{
  Normalize columns of a numeric matrix by robustly estimating a common
  quantile distribution. Unlike standard quantile normalization, columns can
  be assigned weights or excluded from computation of the normalization
  distribution.
}
\usage{
normalize.quantiles.robust(x, copy = TRUE, weights = NULL,
  remove.extreme = c("variance", "mean", "both", "none"),
  n.remove = 1, use.median = FALSE, use.log2 = FALSE,
  keep.names = FALSE)
}
\arguments{
  \item{x}{A numeric matrix, with rows representing probes or features and
    columns representing chips or samples. Missing values are not supported;
    the matrix should contain no \code{NA} or infinite values.}

  \item{copy}{Logical. If \code{TRUE}, normalize a copy of \code{x} and leave
    the input matrix unchanged. If \code{FALSE}, normalize the matrix in place
    to reduce memory use.}

  \item{weights}{Controls construction of the normalization distribution. It
    can be a numeric vector of length \code{ncol(x)} containing non-negative
    weights for the columns, or \code{"huber"} to use Huber M-estimation at
    each quantile with tuning constant 1.345. If \code{NULL}, weights are
    generated from \code{remove.extreme} and \code{n.remove}.}

  \item{remove.extreme}{Character string specifying which columns are assigned
    zero weight when \code{weights} is \code{NULL}. \code{"variance"} removes
    columns with the most extreme variance, \code{"mean"} removes columns whose
    mean is most different from the other columns, \code{"both"} applies both
    checks, and \code{"none"} assigns equal weight to every column. Ignored when
    \code{weights} is supplied.}

  \item{n.remove}{Number of extreme columns to remove for each selected
    criterion. Must be a non-negative integer.}

  \item{use.median}{Logical. If \code{TRUE}, use the median of the sorted
    values at each rank as the normalization distribution. If \code{FALSE},
    construct a weighted mean or Huber estimate, depending on \code{weights}.}

  \item{use.log2}{Logical. If \code{TRUE}, estimate the normalization
    distribution on the log2 scale and transform it back to the original scale.
    This produces a geometric rather than arithmetic mean for the weighted
    normalization. Has no effect when \code{use.median = TRUE}.}

  \item{keep.names}{Logical. If \code{TRUE}, assign the row and column names of
    \code{x} to the normalized matrix.}
}
\details{
  Robust quantile normalization generalizes standard quantile normalization by
  reducing the influence of unusually variable or unusually shifted columns
  when estimating the common reference distribution. For each rank, sorted
  values from the columns are combined according to the requested weighting
  scheme. Every normalized column is then replaced with values from this
  common distribution in rank order.

  By default, \code{remove.extreme} is used to identify columns that should not
  contribute to the reference distribution. With numeric weights, positive
  values give columns greater influence, and zero values exclude columns from
  computation of the reference distribution. Supplying \code{"huber"} uses
  robust Huber estimation at each rank instead of fixed column weights.

  With \code{use.median = TRUE}, the median target is computed independently
  of \code{weights}; this option takes precedence over the weighted mean and
  Huber options. Log2 processing estimates either the weighted mean or Huber
  location on the log2 scale before transforming the target distribution back
  to the intensity scale.

  The implementation does not support missing values. In particular, NA
  handling differs from \code{normalize.quantiles}; do not use this routine
  unless all entries of \code{x} are observed and finite.

  If \code{copy = FALSE}, the input matrix is normalized in place and the
  returned matrix may be the same R object as \code{x}. Use the default
  \code{copy = TRUE} when the original data must be preserved.
}
\value{
  A numeric matrix with the same dimensions as \code{x}. When
  \code{keep.names = TRUE}, row and column names are preserved from \code{x}.
}
\note{
  This function is experimental. The extreme-column selection and robust target
  estimates can be sensitive to small sample sizes.
}
\references{
  Bolstad, B. M., Irizarry, R. A., Astrand, M., and Speed, T. P. (2003).
  A comparison of normalization methods for high density oligonucleotide array
  data based on bias and variance. \emph{Bioinformatics} 19(2), 185--193.
}
\seealso{
  \code{\link{normalize.quantiles}}
}
\examples{
set.seed(1234)

# Create positive intensity-like data with one clearly unusual column.
x <- matrix(rgamma(100, shape = 2, rate = 1), nrow = 20, ncol = 5)
x[, 5] <- x[, 5] * 3

# Assign names before demonstrating the keep.names option.
rownames(x) <- paste0("feature", seq_len(nrow(x)))
colnames(x) <- paste0("sample", seq_len(ncol(x)))

# By default, the function automatically gives zero weight to an extreme
# column selected using the requested remove.extreme criterion.
y <- normalize.quantiles.robust(x, remove.extreme = "variance")
stopifnot(all(dim(y) == dim(x)))

# Supply fixed non-negative column weights. This gives the fifth column only
# half as much influence as the other columns in the reference distribution.
weights <- c(1, 1, 1, 1, 0.5)
y.weights <- normalize.quantiles.robust(x, weights = weights)

# Huber weighting applies robust weights row by row when forming the target
# distribution. The tuning constant is fixed internally at 1.345.
y.huber <- normalize.quantiles.robust(x, weights = "huber")

# Median normalization uses the median value at each rank and does not use
# the fixed weights supplied through the weights argument.
y.median <- normalize.quantiles.robust(x, weights = weights,
                                     use.median = TRUE)

# Fit the robust target distribution on the log2 scale and transform it back
# to the intensity scale.
y.log2 <- normalize.quantiles.robust(x, weights = "huber", use.log2 = TRUE)

# Preserve the input dimnames in the returned matrix.
y.named <- normalize.quantiles.robust(x, keep.names = TRUE)
identical(colnames(y.named), colnames(x))

# copy = FALSE normalizes x in place. This option is useful for memory savings
# when the original matrix is not needed.
x.in.place <- x
y.in.place <- normalize.quantiles.robust(x.in.place, copy = FALSE)
identical(y.in.place, x.in.place)
}
\author{Ben Bolstad, \email{bmb@bmbolstad.com}}
\keyword{manip}




