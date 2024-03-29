% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_fig_one.R
\name{create_fig_one}
\alias{create_fig_one}
\title{Creates Simulation Figure 1}
\usage{
create_fig_one(num_sims = 5, K_f = 2:8, L_f = 4, K = 2:8, L = 4,
  J = c(50, 50, 50), R = 5, f_s = c(function(x) 3 * ((1/2) *
  dbeta(x/R, 1, 8) + (1/2) * dbeta(x/R, 6, 1)), function(x) 3 * ((1/5) *
  dbeta(x/R, 3, 2) + (2/3) * dbeta(x/R, 3, 1) + (2/15) * dbeta(x/R, 1, 1)),
  function(x) 3 * ((1/2) * dbeta(x/R, 8, 2) + (1/2) * dbeta(x/R, 30, 50))),
  iter_max = 15000, warm_up = 12000, thin = 1, seed = NULL)
}
\arguments{
\item{K_f}{*vector* of the different number of cluster components in fhmm}

\item{L_f}{number of within-cluster mixture components in fhmm}

\item{K}{number of cluster sticks for NDP}

\item{L}{number of within-cluster mixture-sticks for NDP}

\item{J}{total number of groups}

\item{f_s}{list of functions that should be the same length as K_true}

\item{iter_max}{maximum number of iterations to use in model fitting}

\item{warm_up}{number of warm-up iterations to use in model fitting}
}
\value{
plot showing results
}
\description{
This function creates a table comparing the novel NDP
to a classic heirarchical finite mixture model
}
