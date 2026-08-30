\name{rcModelPLMd}
\alias{rcModelPLMd}

\title{Fit row-column models with group-specific row effects}
\description{
  These functions fit additive row-column decomposition models in which the
  row effects are allowed to vary by predefined groups of arrays. This is useful
  when a subset of probes behaves differently across blocks, batches, or other
  sample groups, and a common row effect is not adequate for all arrays. The
  procedure automatically checks whether selected rows show evidence of group-
  dependent behavior and, if so, splits the row effect into separate values for
  each group.
}
\usage{
rcModelPLMd(y, group.labels)
}
\arguments{
  \item{y}{A numeric matrix with rows representing probes/features and columns
    representing arrays/samples. The model is
    \deqn{y_{ij} = r_i + c_j + \epsilon_{ij}}{y_ij = r_i + c_j + e_ij},
    where \eqn{r_i} is the row effect and \eqn{c_j} is the column effect.}

  \item{group.labels}{A vector or factor specifying the group label for each
    column of \code{y}. Its length must be \code{ncol(y)}. Each group must
    contain at least two arrays. The procedure tests whether row effects differ
    across these groups and, where appropriate, splits a row effect into group-
    specific effects.}
}
\value{
  A list with the following components:
  \describe{
    \item{Estimates}{The fitted parameter estimates. These include the column
      effects and the row effects, with additional group-specific row effects
      for any rows whose effects are split across the groups. The exact order is
      the same as the fitted design matrix: column effects first, followed by the
      row effects and any split terms needed for the group-specific rows.}

    \item{Weights}{The final observation weights used during the robust fit.}

    \item{Residuals}{A matrix of residuals after subtracting the fitted row and
      column effects.}

    \item{StdErrors}{The standard errors associated with the fitted parameters,
      in the same order as \code{Estimates}.}

    \item{WasSplit}{An integer vector of length \code{nrow(y)} indicating which
      rows were split into separate effects across groups. A value of 1 means
      that the row effect was allowed to differ across the specified groups; 0
      means it remained shared across all groups.}
  }
}
\details{
  The ordinary PLM model fits a single shared row effect for all arrays in the
  matrix, as in \deqn{y_{ij} = r_i + c_j + \epsilon_{ij}}{y_ij = r_i + c_j + e_ij}.
  The PLM-d procedure generalizes this model by allowing some row effects to
  vary across group-defined blocks. Thus, for a row that appears to behave
  differently in one group relative to another, the model estimates separate row
  effects for each group instead of a single common value.

  The implementation in the C code follows the idea described in the PLM
  extension manuscript: the method first fits a standard robust row-column model,
  then examines the residuals for each row to see whether they are associated
  with the grouping variable. If a row shows evidence of a group-specific
  departure, that row is split and fitted with separate row effects for the
  different groups. The process repeats until no further group-dependent row
  effects are detected.

  This differs from the robust row/column weighting methods in
  \code{\link{rcModelPLMr}}. The PLM-r family addresses outlying rows or columns
  by downweighting them during the fit. The PLM-d family instead asks whether the
  underlying row effect itself differs across groups and, if so, estimates a
  separate effect within each group.
}
\section{Model Overview}{
  The PLM-d model retains the additive form
  \deqn{y_{ij} = r_i + c_j + \epsilon_{ij}}{y_ij = r_i + c_j + e_ij}
  but allows selected rows to be split into group-specific effects. In practice,
  this is useful when a subset of probes has a different systematic offset in one
  biological or experimental group than in another. Rows that do not show such a
  difference keep a single shared effect across all groups.
}
\section{Choosing a Method}{
  \describe{
    \item{\code{\link{rcModelPLM}}}{Use this for a standard robust decomposition
      with one shared row effect across all arrays. This is the simplest baseline
      model when no grouping structure is expected.}

    \item{\code{\link{rcModelPLMr}}}{Use this when the main issue is rows or
      columns that are unusually noisy or outlying. These methods downweight such
      rows or columns rather than splitting the row effect by group.}

    \item{\code{rcModelPLMd}}{Use this when you have a known grouping variable and
      suspect that only some rows behave differently between groups. This method
      tests for and fits group-specific row effects only where the data support
      them.}
  }
}
\section{Comparison to ordinary PLM and PLM-r models}{
  \describe{
    \item{\code{\link{rcModelPLM}}}{Fits a common row effect across all columns.
      There is no notion of group-specific row effects.}

    \item{\code{\link{rcModelPLMr}}}{Adds row and/or column robustness by
      downweighting outlying rows or columns. It does not split row effects by
      group labels.}

    \item{\code{rcModelPLMd}}{Retains the additive structure but allows selected
      probes to have separate effects in different groups when the data support
      that structure. This is a group-wise model extension rather than a
      weighting-based robustification.}
  }
}
\seealso{\code{\link{rcModelPLM}},\code{\link{rcModelPLMr}},\code{\link{rcModelPLMrr}},\code{\link{rcModelPLMrc}}}
\references{
  Bolstad, B. M. PLM Extensions Manuscript.
  \url{https://bmbolstad.com/stuff/PLM-Extensions-Manuscript.pdf}
}
\examples{
set.seed(1234)

# Example 1: two groups with a row-level effect that differs between them.
# This matches the intended PLM-d use case: the algorithm should detect a
# group-specific row effect and split it.
row.effects <- c(4, 3, 2, 1, -1, -2, -3, -4)
col.effects <- c(8, 9, 10, 11, 12, 10)
y <- outer(row.effects, col.effects, "+") + rnorm(48, sd = 0.25)

y[8, 4:6] <- c(11, 12, 10) + 2.5 + rnorm(3, sd = 0.25)
y[5, 4:6] <- c(11, 12, 10) - 2.5 + rnorm(3, sd = 0.25)

fit <- rcModelPLMd(y, group.labels = factor(c(1, 1, 1, 2, 2, 2)))
fit$WasSplit
fit$Estimates

# Example 2: the same design without any group-specific row shift.
# This checks the non-split case, where all row effects remain common across groups.
row.effects <- c(seq(-0.5, -0.1, by = 0.1), seq(0.1, 0.5, by = 0.1))
y.clean <- outer(row.effects, col.effects, "+") + rnorm(60, sd = 0.1)
fit.clean <- rcModelPLMd(y.clean, group.labels = factor(c(1, 1, 1, 2, 2, 2)))
fit.clean$WasSplit

# Compare the residual structure with ordinary PLM.
par(mfrow = c(2, 2))
matplot(y, type = "l", col = c(rep("red", 3), rep("blue", 3)),
        ylab = "signal", xlab = "probe", main = "Observed data")
matplot(rcModelPLM(y)$Residuals, col = c(rep("red", 3), rep("blue", 3)),
        ylab = "residuals", xlab = "probe", main = "Residuals (PLM)")
matplot(fit$Residuals, col = c(rep("red", 3), rep("blue", 3)),
        ylab = "residuals", xlab = "probe", main = "Residuals (PLM-d)")
}
\author{B. M. Bolstad \email{bmb@bmbolstad.com}}
\keyword{models}
