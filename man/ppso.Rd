\name{ppso}
\alias{ppso}
\title{ppso-package, general information}
\description{
  General information about the ppso package.
}
\details{
  (Optionally parallelized) optimization using PSO (Particle Swarm Optimzation) or DDS (Dynamically Dimensioned Search) algorithms, which excell for multidimensional (multi-parameter) functions with many local extrema and a restricted number of function evaluations. The parallelized version builds on Rmpi and is intended for highly computationally intensive objective functions (>20 s evaluation time). Another focus of this package is the possibility to resume interrupted optimization runs from the intermediate project files. These features make this package useful for the automatic calibration of complex numerical models (e.g. hydrological models).
-- experimental version, use with care. Comments welcome --
  Type \code{library(help = ppso)} to see a list of (public) ppso functions.
}
\author{Till Francke}
