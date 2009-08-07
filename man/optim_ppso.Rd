\name{optim_pso}
\Rdversion{1.1}
\alias{optim_pso}
\alias{optim_ppso}
\alias{optim_ppso_robust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
optim_pso
}
\description{
Optimize a given function \code{objective_function} with regard to its parameters contained in the vector \code{par} towards a minimum value of f
}
\usage{
optim_pso        (objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4), load_projectfile="try", break_file=NULL)
optim_ppso       (objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4), load_projectfile="try", break_file=NULL, nslaves = 3)
optim_ppso_robust(objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4), load_projectfile="try", break_file=NULL, nslaves = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{


  \item{objective_function}{
function whose numerical return value is to be minimized.  \code{objective_function} needs to accept a vector \code{X} as argument, which contains \code{number_of_parameters} numerical values.
\code{objective_function} is expected to return a numerical value and handle exceptions (out-of-bounds-parameters, etc.) on its own by returning \code{Inf}.
}
  \item{number_of_parameters}{
Number of parameters to be optimized.
}
  \item{number_of_particles}{
Number of particles used in the Particle Swarm Optimization.
}
  \item{max_number_of_iterations}{
Abortion criteria: Maximum number of iterations (i.e. update of particle positions) that are used.
}
  \item{w}{
Inertia constant, i.e. "weight" of the particles.
}
  \item{C1}{
Cognitive component. Weighting factor for the "personal" experience of each particle.
}
  \item{C2}{
Social component. Weighting factor for the "swarm" experience of each particle.
}
  \item{abstol}{
Abortion criteria: minimum absolute improvement between iterations  (default: -Inf)
}
  \item{reltol}{
Abortion criteria: minimum absolute improvement between iterations (default: -Inf)
}
  \item{max_wait_iterations}{
Number of iterations, within these an improvement of the above described (\code{abstol, reltol}) quality has to be achieved before abortion
}
  \item{wait_complete_iteration}{
TRUE:  Update particle velocities after function \code{objective_function} has been evaluated for ALL particles.
FALSE: Update particle velocities after each completed evaluation of function \code{objective_function} (suggested in parallel mode).
}
  \item{parameter_bounds}{
matrix containing lower (first column) and upper (second column) boundary for each parameter
}
  \item{Vmax}{
maximum velocity of particles in parameter space
}

  \item{lhc_init}{
set starting positions of particles using Latin Hypecube Sampling (TRUE, requires package \code{lhc}) or purely random (FALSE)
}

  \item{logfile}{
Name of logfile for optional logging of ALL model runs. NULL disables logging. If an prior optimization run is successfully resumed (see \code{projectfile}), an existing \code{logfile} is appended to.
}

  \item{projectfile}{
Name of project file for optional logging of state of optimization, which enables resuming aborted optimizations. An existing file is overwritten! NULL disables project file.
}

  \item{save_interval}{
minimum number of function evaluations to compute before the \code{projectfile} is updated.
}

  \item{load_projectfile}{
 \code{"yes"}: require \code{projectfile} to initialise particles;  \code{"try"}: load \code{projectfile} , if existent;  \code{"no"}: ignore existing \code{projectfile}
}

  \item{break_file}{
Name of a file file that will cause the termination of the optimization when it is encountered. Useful for interrupting a running optimization gracefully. At start-up, a existing file of this name will be deleted.
}


for optim_ppso*:
  \item{nslaves}{
number of rmpi slaves to spawn (default -1: as many as possible, requires package \code{Rmpi})
}

(didactical purpose, only for two-parameter search and fast objective function. "base" works only with  wait_complete_iteration=TRUE)
  \item{do_plot}{
enable 3D-plot of response surface and search progress.
"base": use basic 3D-plotting functions.
"rgl":  use moveable rgl-plotting (requires package \code{rgl})
}
  \item{wait_for_keystroke}{
waiting for keystroke between iterations (e.g. for plotting)
}
}


\details{
recommended use:
\code{optim_pso}: serial computation for fast objective functions and single core machines
\code{optim_ppso}			: parallel computation: computationally demanding objective function and multi-core CPU or network slaves.       
\code{optim_ppso_robust}	: slightly slower than \code{optim_ppso}, but presumably more robust (to be tested), since it doesn't leave slaves in "listen" mode after crashes

}
\value{
\code{ret_val=list(value,par,iterations)}
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{value }{value of objective function at minimum}
  \item{par }{parameters of minimum found}
  \item{iterations }{number of iterations needed}
  \item{break_flag }{criterium that caused the termination of the algorithm (\code{c("max iterations reached", "converged", "all slaves gone", "user interrupt")})}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Till Francke <francke_at_uni-potsdam.de>
}
\note{
- still in development, comments welcome -
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

optim_pso(objective_function=sample_function)



set.seed(2)
optim_ppso_robust(objective_function=sample_function,projectfile=NULL)

set.seed(2)
optim_ppso_robust(objective_function=sample_function,projectfile=NULL,nslaves=1)


set.seed(2)
optim_ppso(objective_function=sample_function)

set.seed(2)
optim_ppso(objective_function=sample_function, abstol=0.01,max_wait_iterations=2,max_number_of_iterations=10,load_projectfile="no")

set.seed(2)
optim_ppso(objective_function=sample_function, abstol=0.1,max_wait_iterations=2,max_number_of_iterations=10,load_projectfile="yes")


optim_ppso(objective_function=sample_function)
optim_ppso_robust(objective_function=sample_function)
}

\keyword{ ~optimization }
\keyword{ ~parallel } 
\keyword{ ~swarm } 
\keyword{ ~particle } % __ONLY ONE__ keyword per line