\name{mix}
\alias{mix}
\alias{Bmix-package}
\alias{Bmix}

\title{ Bayesian inference for (dynamic) stick-breaking mixtures }

\description{ Basic implementation for a variety of stick-breaking
mixture models, including particle learning and Gibbs sampling for
static DP mixtures, particle learning for dynamic BAR stick-breaking,
and DP mixture regression.  }

\details{

This is a bare-bones implementation of sampling algorithms for
Bayesian stick-breaking mixture models.  The software is designed to
be easy to customize to suit different situations and for
experimentation with stick-breaking models.  Since particles are
repeatedly copied, it is not especially efficient.

The package implements particle learning (Carvalho et al, 2009) for
both dynamic and constant stick-breaking mixture models, and collapsed
Gibbs sampling for DP mixtures.  Conditional sufficient statistics for
each mixture component are output as \sQuote{particle} files.

Mixture kernels are the product of independent multinomial densities
for each categorical variable, and a multivariate normal density for
continuous covariates.  The base measure is conditionally conjugate
normal-Wishart-Dirichlet product, with Wishart hyperprior for inverse base
covariance.  Beta-autoregressive stick-breaking is used to model
correlated densities.

Refer to Taddy (2009) for all specification details.  
 
See DPreg demo for regression with categorical and continuous covariates, 
with additional Gibbs sampling for filtered particles. 

See bar1D and bar2D demos for dynamic stick-breaking mixture density estimation, 
Bayes factor calculations, and comparison between correlated and independent model fit.  
}

\usage{
mix(Z, alpha, g0params, times=NULL,  rho=NULL,  cat=0, 
    state=NULL, read=FALSE, print=FALSE, N=100, niter=0) }

\author{
Matt Taddy:
\email{matt.taddy@chicagobooth.edu}
}
\references{ 

\emph{ An auto-regressive mixture model for dynamic spatial Poisson processes: 
Application to tracking the intensity of violent crime} (Taddy 2009), 

\emph{Particle learning for general mixtures} (Carvalho, Lopes, Polson, and Taddy 2009), 

\emph{ A Bayesian nonparametric approach to inference for quantile regression} (Taddy and Kottas 2009). 

and other papers at \url{faculty.chicagobooth.edu/matt.taddy/research.}}


\arguments{
  \item{Z}{\code{data.frame} of observations, with the last \code{cat} columns categorical variables. }
  \item{alpha}{ Stick-breaking precision parameter. }
  \item{g0params}{ A vector of the base measure parameters:
   \eqn{[\gamma, \kappa, \nu, \gamma_\Omega, \Psi_\Omega ]},
  as in Taddy (2009), followed by a list of Dirichlet parameters for each categorical variable.}
  \item{times}{ Discrete time period for each data row; must be an increasing list of integers starting at 0.}
  \item{rho}{ Stick-breaking correlation parameter (between 0 and 1). }
  \item{cat}{ Number of categorical variables.}
  \item{state}{ Random number generator seed. }
  \item{read}{ Read in existing particle info? }
  \item{print}{ Print out particles for each new time? WARNING: This will litter your working directory with .particle* files!}
  \item{N}{ Number of particles.}
  \item{niter}{ Number of Gibbs sampling iterations after filtering (only for \code{times = NULL}). }
}

\value{
 
If print=TRUE particle representations for each time point are printed
to files .particle[i].[t].[rho].txt, where \sQuote{i} is the particle
id, \sQuote{t} is the time, and \sQuote{rho} is the stick-breaking correlation
parameter.  These files are read into R via the \code{particle()} function.

Output to R is minimal.  The object returned by function \code{mix()} 
is mostly just a list of input variables, except for

\item{logprob}{Filtered marginal log-likelihood}
\item{m}{Filtered mean number of mixture components}
\item{k}{If kout, this is the list of filtered allocations}

}

\seealso{\code{particle}}

\examples{
demo(package="Bmix")
  }

