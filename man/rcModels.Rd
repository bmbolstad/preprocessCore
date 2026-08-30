\name{rcModels}
\alias{rcModelPLM}
\alias{rcModelWPLM}
\alias{rcModelMedianPolish}
\title{Fit additive row-column decomposition models}
\description{Decompose a matrix into additive row and column effects plus residuals. 
These functions fit row-column models useful for extracting systematic patterns and 
identifying systematic biases in data (e.g., probe effects and sample effects in microarray data). 
Three fitting methods are provided: robust linear models with and without observation weights, 
and median polish for a non-parametric alternative.
}
\usage{
rcModelPLM(y,row.effects=NULL,input.scale=NULL)
rcModelWPLM(y, w,row.effects=NULL,input.scale=NULL)
rcModelMedianPolish(y)
}
\arguments{
  \item{y}{A numeric matrix with rows representing probes/features and columns 
    representing arrays/samples. Used by all three functions to fit the row-column model.}
  \item{w}{(Only for rcModelWPLM) A matrix or vector of non-negative weights for each observation.
    If a vector, must have length \code{prod(dim(y))} (i.e., total number of elements).
    If a matrix, must have the same dimensions as \code{y}. Zero weights completely exclude 
    observations from the fit; larger weights increase an observation's influence on estimation.}
  \item{row.effects}{(Only for rcModelPLM and rcModelWPLM) A vector of pre-computed row effects 
    to use instead of fitting them. Must have length \code{nrow(y)} and must sum to exactly zero 
    (within machine precision). When supplied, the function fixes these row effects and estimates 
    only column effects, analyzing each column separately. Useful for applying probe effects from 
    a reference dataset to new data or for sensitivity analysis. If \code{NULL} (default), row 
    effects are estimated from the data.}
  \item{input.scale}{(Only for rcModelPLM and rcModelWPLM) A pre-computed scale (residual standard 
    deviation) estimate to use instead of estimating it from the data. Must be positive. Useful 
    for consistency across multiple analyses or when external scale estimates are available. 
    If \code{NULL} (default), scale is estimated from the residuals using a robust method.}
}
\value{
  A list with the following components. Components are named and can be accessed via 
  \code{result$Estimates}, \code{result$Residuals}, etc.
  
  \describe{
    \item{Estimates}{A numeric vector of length \code{nrow(y) + ncol(y)} containing 
      the fitted row and column effects. Elements are ordered as: 
      [column effect 1, ..., column effect p, row effect 1, ..., row effect n] 
      where n = \code{nrow(y)} and p = \code{ncol(y)}. These are the constrained 
      parameters (sum to zero for rcModelPLM/rcModelWPLM; median = 0 for rcModelMedianPolish).}
    
    \item{Weights}{For rcModelPLM and rcModelWPLM: A numeric matrix (n × p) of final weights 
      used in fitting, where weights < 1 indicate downweighting of outliers. 
      For rcModelMedianPolish: \code{NULL} (weights not computed for median polish).}
    
    \item{Residuals}{A numeric matrix (n × p) of residuals after subtracting fitted effects 
      from the data. For a good fit, residuals should appear random with no systematic patterns.}
    
    \item{StdErrors}{For rcModelPLM and rcModelWPLM: A numeric vector of length 
      \code{nrow(y) + ncol(y)} containing standard errors for each parameter estimate, 
      in the same order as Estimates. Larger values indicate less precise estimation. 
      For rcModelMedianPolish: \code{NULL} (standard errors not computed; median polish 
      is non-parametric and does not provide SEs).}
    
    \item{Scale}{For rcModelPLM and rcModelWPLM: A numeric scalar (length 1) containing 
      the estimated scale parameter (residual standard deviation). Indicates the overall 
      noise level in the data. For rcModelMedianPolish: This component is not returned 
      (median polish does not estimate scale).}
  }
}
\details{
  These functions fit row-column models to the specified input
  matrix. Specifically the model \deqn{y_{ij} = r_i + c_j +
    \epsilon_{ij}}{y_ij = r_i + c_j + e_ij}
       with \eqn{r_i} and \eqn{c_j} as row and column effects
       respectively. Note that this functions treat the row effect as
       the parameter to be constrained using sum to zero (for
       \code{rcModelPLM} and \code{rcModelWPLM}) or median of zero (for
       \code{rcModelMedianPolish}).

       The \code{rcModelPLM} and \code{rcModelWPLM} functions use a
       robust linear model procedure for fitting the model.

       The function \code{rcModelMedianPolish} uses the median polish algorithm.
}
\section{Model Overview}{

  These functions decompose a matrix into additive row and column effects plus residuals 
  using the linear model: \deqn{y_{ij} = r_i + c_j + \epsilon_{ij}}{y_ij = r_i + c_j + e_ij}
  
  This decomposition is useful for:
  \describe{
    \item{Detecting systematic patterns}{Identify whether variation in your data 
      is systematic (captured by row/column effects) or random (residuals).}
    
    \item{Quality assessment}{Row effects represent probe/feature-specific biases; 
      column effects represent array/sample-specific effects. Unusual values can 
      indicate quality issues.}
    
    \item{Normalization and preprocessing}{Understanding and removing systematic effects 
      before downstream analysis.}
    
    \item{Comparative analysis}{Compare the magnitudes of row vs. column effects to 
      understand whether variation is more probe-specific or sample-specific.}
  }
}
\section{Choosing a Method}{

  \describe{
    \item{rcModelPLM (Robust Linear Model)}{
      Use this for standard row-column decomposition. It employs Huber's robust 
      M-estimation procedure (tuning constant k=1.345) that downweights 
      outliers. Good for data with occasional extreme values. Estimates row effects 
      that sum to exactly zero (identifiability constraint). Returns scale estimate and 
      standard errors.
    }
    
    \item{rcModelWPLM (Weighted Robust Linear Model)}{
      Use this when you have prior weights for each observation (e.g., quality weights, 
      known measurement uncertainties). Weights must be non-negative; zero weights 
      exclude observations entirely. Uses the same Huber's M-estimation procedure as rcModelPLM 
      but incorporates the weights. Useful for combining data from multiple sources with 
      different reliabilities.
    }
    
    \item{rcModelMedianPolish}{
      Use this for a non-parametric alternative to RLM. More resistant to outliers than 
      even robust RLM, and makes fewer distributional assumptions. Iteratively refines 
      row and column medians. Constraints row effects to have median equal to zero. 
      Does not return scale estimate or standard errors. Simpler interpretation but 
      less statistical efficiency than RLM when data is approximately normal.
    }
  }
}
\section{Using Pre-Computed Row Effects}{

  You can provide pre-computed row effects (e.g., from a reference dataset) using the 
  \code{row.effects} parameter. When supplied:
  \describe{
    \item{Constraints}{Row effects must sum to exactly zero (within machine precision). 
      If this constraint is violated, the function will error.}
    
    \item{Behavior}{The model fits column effects given the fixed row effects, and 
      analyzes each column separately. This is useful for ensuring consistency across 
      multiple analyses.}
    
    \item{Use cases}{Applying probe effects from a reference or training dataset to new data; 
      sensitivity analysis; two-stage model fitting.}
  }
}
\section{Interpreting the Output}{

  The output is a list with components:
  \describe{
    \item{Estimates}{Concatenated vector of column effects followed by row effects. 
      This order is chosen for historical reasons and model fitting efficiency.}
    
    \item{StdErrors}{Standard errors for each parameter in the same order as Estimates. 
      Indicates precision of each effect estimate. Only returned by rcModelPLM and 
      rcModelWPLM (rcModelMedianPolish returns NULL).}
    
    \item{Scale}{Estimated scale parameter (residual variability). Only returned by 
      rcModelPLM and rcModelWPLM. Useful for understanding noise level in data.}
    
    \item{Residuals}{The residuals \eqn{e_ij}{e_ij} after subtracting estimated effects. 
      Should appear random if model fits well.}
    
    \item{Weights}{The final observation weights used in fitting. For rcModelPLM, 
      these are computed based on residuals; for rcModelWPLM, combines input weights 
      with residual-based weights.}
  }
}
\section{Related Functions}{
  \code{\link{rcModelPLMr}} and \code{\link{rcModelPLMd}} are related variants that fit 
  row-column models with additional row-specific constraints. These functions constrain 
  row effects to be equal across all columns (rcModelPLMr) or use row-specific digestion 
  effects (rcModelPLMd), respectively. See their documentation for details on when to use 
  these variants versus the standard rcModel functions.
}
\seealso{\code{\link{rcModelPLMr}},\code{\link{rcModelPLMd}}}
\examples{
# Set random seed for reproducible examples
set.seed(1234)

# Create sample data: 10 row effects (probes) and 5 column effects (arrays)
row.effects <- c(seq(-0.5, -0.1, by = 0.1), seq(0.1, 0.5, by = 0.1))
col.effects <- c(10, 11, 10.5, 12, 9.5)

# Perfect model: no noise
y.perfect <- outer(row.effects, col.effects, "+")

# Noisy data: add random noise to model
y <- y.perfect + rnorm(50)

# Create weights for weighted model (optional)
w <- runif(50)

# Basic usage: fit all three methods without additional parameters
rcModelPLM(y)
rcModelWPLM(y, w)
rcModelMedianPolish(y)

# Using pre-computed row effects (e.g., from a reference dataset)
# When row effects are supplied, only column effects are estimated
rcModelPLM(y, row.effects = row.effects)
rcModelWPLM(y, w, row.effects = row.effects)

# Using a fixed input scale parameter
rcModelPLM(y, input.scale = 1.0)
rcModelWPLM(y, w, input.scale = 1.0)

# Combining multiple options: fixed row effects AND fixed scale
rcModelPLM(y, row.effects = row.effects, input.scale = 1.0)
rcModelWPLM(y, w, row.effects = row.effects, input.scale = 1.0)

# Examine the output structure
result <- rcModelPLM(y)
str(result)  # View the structure
result$Estimates  # Extract column and row effects (column effects first)
dim(result$Residuals)  # Residuals are same dimension as input
}
\author{B. M. Bolstad \email{bmb@bmbolstad.com}}
\keyword{models}
