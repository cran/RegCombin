% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_stat_variance.R
\name{compute_stat_variance}
\alias{compute_stat_variance}
\title{Function to compute the Variance bounds on the noncommon regressor Xnc}
\usage{
compute_stat_variance(
  sample1 = NULL,
  X1_x,
  X2,
  X1_y,
  Y,
  values,
  refs0,
  dimX1,
  dimX2,
  nb_pts,
  sam0,
  lim = 1,
  weights_x = NULL,
  weights_y = NULL,
  constraint = NULL,
  c_sign = NULL,
  nc_sign = NULL,
  values_sel = NULL
)
}
\arguments{
\item{sample1}{if NULL compute the point estimate, if a natural number then evaluate a bootstrap or subsampling replication.}

\item{X1_x}{the common regressor on the dataset  (Xnc,Xc). Default is NULL.}

\item{X2}{the noncommon regressor on the dataset  (Xnc,Xc). No default.}

\item{X1_y}{the common regressor on the dataset  (Y,Xc). Default is NULL.}

\item{Y}{the outcome variable. No default.}

\item{values}{the different unique points of support of the common regressor Xc.}

\item{refs0}{indicating the positions in the vector values corresponding to the components of betac.}

\item{dimX1}{the dimension of the common regressors Xc.}

\item{dimX2}{the dimension of the noncommon regressors Xnc.}

\item{nb_pts}{the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.}

\item{sam0}{the directions q to compute the variance bounds on the radial function.}

\item{lim}{the limit number of observations under which we do no compute the conditional variance.}

\item{weights_x}{the sampling weights for the dataset (Xnc,Xc).}

\item{weights_y}{the sampling weights for the dataset (Y,Xc).}

\item{constraint}{a vector indicating the different constraints in a vector of the size of X_c indicating the type of constraints, if any on f(X_c) : "concave", "concave", "nondecreasing", "nonincreasing", "nondecreasing_convex", "nondecreasing_concave",  "nonincreasing_convex", "nonincreasing_concave", or NULL for none. Default is NULL, no contraints at all.}

\item{c_sign}{sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.}

\item{nc_sign}{sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.}

\item{values_sel}{the selected values of Xc for the conditioning. Default is NULL.}
}
\value{
a list containing:

 - upper: the upper bound in the specified directions, possibly with sign constraints

 - lower: the lower bound in the specified directions, possibly with sign constraints

 - unconstr: the bounds without sign constraints in the specified directions

 - Ykmean: the means of Y|Xc for the considered sample

 - Xkmean: the means of Xnc|Xc for the considered sample

 - DYk: the difference of means of Y|Xc =k -  Y|Xc =0 for the considered sample

 - DXk: the difference of means of Xnc|Xc =k -  Xnc|Xc =0 for the considered sample

 - tests: the pvalues of the tests H0 : DXk =0

 - ratio_ref: the ratio R in the radial function computed for the initial sample
}
\description{
Function to compute the Variance bounds on the noncommon regressor Xnc
}
