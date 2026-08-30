\name{rcModelPLMr}
\alias{rcModelPLMr}
\alias{rcModelPLMrr}
\alias{rcModelPLMrc}
\alias{rcModelWPLMr}
\alias{rcModelWPLMrr}
\alias{rcModelWPLMrc}

\title{Fit row-column models with row and/or column robustness}
\description{
  These functions fit additive row-column decomposition models while adding an
  additional layer of robustness beyond the standard \code{\link{rcModelPLM}}
  procedure. They are useful when systematic problems affect entire rows or
  columns, rather than only isolated elements in the matrix. The methods differ
  in which direction of robustification is used: both row and column effects,
  row effects only, or column effects only.
}
\usage{
rcModelPLMr(y)
rcModelPLMrr(y)
rcModelPLMrc(y)
rcModelWPLMr(y, w)
rcModelWPLMrr(y, w)
rcModelWPLMrc(y, w)
}
\arguments{
  \item{y}{A numeric matrix with rows representing features or probes and columns
    representing arrays or samples. The model is \\deqn{y_{ij} = r_i + c_j +
    \\epsilon_{ij}}{y_ij = r_i + c_j + e_ij}, where \eqn{r_i} denotes the
    row effect and \eqn{c_j} the column effect.}

  \item{w}{For the weighted variants \code{rcModelWPLMr}, \code{rcModelWPLMrr},
    and \code{rcModelWPLMrc}, a matrix or vector of non-negative observation
    weights. If a vector is supplied, it must have length
    \code{prod(dim(y))}; if a matrix is supplied, it must have the same
    dimensions as \code{y}. Zero weights remove an observation from the fit,
    while larger weights give that observation more influence.}
}
\value{
  A list with the following components:
  \describe{
    \item{Estimates}{A numeric vector of length \code{nrow(y) + ncol(y)}
      containing the fitted row and column effects. The values are stored in
      column-effect then row-effect order: [column effect 1, ..., column effect
      p, row effect 1, ..., row effect n].}

    \item{Weights}{A numeric matrix with the same dimensions as \code{y}. Its
      entries are the final observation weights used in the fit; values near 0
      indicate downweighting of outlying observations.}

    \item{Residuals}{A numeric matrix of residuals after subtracting the fitted
      additive row and column effects.}

    \item{StdErrors}{A numeric vector of standard errors for the fitted effects,
      ordered in the same way as \code{Estimates}.}
  }
}
\details{
  These functions fit the row-column model
  \deqn{y_{ij} = r_i + c_j + \epsilon_{ij}}{y_ij = r_i + c_j + e_ij}
  with row and column effects constrained so that the row effects are treated as
  the identifiable component. The robust fitting is based on the same Huber
  M-estimation framework as \code{\link{rcModelPLM}}, with tuning constant
  \eqn{k = 1.345}, but in addition the methods estimate row and/or column
  robustness weights to downweight systematic row-level or column-level
  departures from the additive structure.

  The three variants differ in how this extra robustness is applied:

  \describe{
    \item{\code{rcModelPLMr} / \code{rcModelWPLMr}}{Fit the PLM-r procedure,
      which adds both row and column robustness to the standard element-wise
      robust fit. This is appropriate when outlying behavior may occur in
      either direction, affecting entire rows or columns of the matrix.}

    \item{\code{rcModelPLMrr} / \code{rcModelWPLMrr}}{Fit the PLM-rr procedure,
      which adds row robustness only. This is useful when a small number of rows
      are systematically discrepant but the columns are not expected to require
      special handling.}

    \item{\code{rcModelPLMrc} / \code{rcModelWPLMrc}}{Fit the PLM-rc procedure,
      which adds column robustness only. This is useful when a small number of
      columns are systematically discrepant but row-specific problems are not the
      main concern.}
  }

  In practice, these methods are most helpful when the data contain structured
  outliers or patterns of contamination that cannot be well represented by the
  ordinary element-wise robust fit alone. The robust weights are designed to
  identify and reduce the influence of rows or columns that deviate from the
  common additive pattern.
}
\section{Comparison to the ordinary PLM functions}{
  The functions \code{\link{rcModelPLM}} and \code{\link{rcModelWPLM}} use a
  robust linear model fit with a single robust weighting step at the element
  level. The PLM-r family extends this idea by adding an additional row/column
  weighting layer. As a result, the PLM-r variants are generally more appropriate
  when the problem is not just a few isolated outlying cells but broader
  row-wise or column-wise departures from the additive form.
}
\section{Comparison to \code{rcModelPLMd}}{
  The function \code{\link{rcModelPLMd}} handles a different kind of structure:
  it tests whether row effects differ across groups and allows a row effect to be
  split by group labels. In contrast, the PLM-r variants do not split effects by
  group labels; instead they identify rows or columns whose residual behavior is
  inconsistent with the overall additive pattern and downweight them during the
  fit.
}
\seealso{\code{\link{rcModelPLM}},\code{\link{rcModelWPLM}},\code{\link{rcModelPLMd}}}
\references{
  Bolstad, B. M. PLM Extensions Manuscript.
  \url{https://bmbolstad.com/stuff/PLM-Extensions-Manuscript.pdf}
}
\examples{
set.seed(1234)

# Create a clean matrix with additive row and column effects.
col.effects <- c(10, 11, 10.5, 12, 9.5)
row.effects <- c(seq(-0.5, -0.1, by = 0.1), seq(0.1, 0.5, by = 0.1))
y <- outer(row.effects, col.effects, "+") + rnorm(50, sd = 0.1)

# The usual row-column fit with row/column robustness.
fit <- rcModelPLMr(y)
fit$Estimates
fit$Weights[1:3, 1:3]

# A row outlier can be handled especially well by the row-robust fit.
y.row <- y
y.row[1, ] <- y.row[1, ] + 3
fit.r <- rcModelPLMr(y.row)
fit.rr <- rcModelPLMrr(y.row)
fit.rc <- rcModelPLMrc(y.row)

# The row-robust variant should reduce the impact of the anomalous row.
summary(c(fit.r$Weights))
summary(c(fit.rr$Weights))
summary(c(fit.rc$Weights))

# A column outlier can be handled especially well by the column-robust fit.
y.col <- y
y.col[, 4] <- y.col[, 4] + 3
fit.r.col <- rcModelPLMr(y.col)
fit.rr.col <- rcModelPLMrr(y.col)
fit.rc.col <- rcModelPLMrc(y.col)

# The column-robust variant should reduce the impact of the anomalous column.
summary(c(fit.r.col$Weights))
summary(c(fit.rr.col$Weights))
summary(c(fit.rc.col$Weights))

# Visualize the final weights for the row-outlier case.
par(mfrow = c(2, 2))
image(1:nrow(y.row), 1:ncol(y.row), fit.r$Weights,
      xlab = "row", ylab = "column", main = "PLM-r", zlim = c(0, 1))
image(1:nrow(y.row), 1:ncol(y.row), fit.rr$Weights,
      xlab = "row", ylab = "column", main = "PLM-rr", zlim = c(0, 1))
image(1:nrow(y.row), 1:ncol(y.row), fit.rc$Weights,
      xlab = "row", ylab = "column", main = "PLM-rc", zlim = c(0, 1))
matplot(y.row, type = "l", main = "Data with a row outlier")

# Weighted variants allow prior knowledge to be incorporated.
w <- rep(1, length(y))
fit.w <- rcModelWPLMr(y, w)
fit.w$Estimates
}
\author{B. M. Bolstad \email{bmb@bmbolstad.com}}
\keyword{models}
