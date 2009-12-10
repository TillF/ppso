\name{optim_pso}
\Rdversion{1.1}
\alias{optim_pso}
\alias{optim_ppso}
\alias{optim_ppso_robust}
\alias{optim_dds}
\alias{optim_pdds_robust}

\title{
(Parallel) optimization using Particle Swarm Optimization or Dynamically Dimensioned Search
}
\description{
\code{optim_pso} optimizes a given function \code{objective_function} with regard to its parameters contained in the vector \code{par} towards a minimum value of f using Particle Swarm Optimization.
\code{optim_ppso} and \code{optim_ppso_robust} are parallelized versions (using multiple CPUs).
\code{optim_dds} optimizes using a Dynamically Dimensioned Search with  \code{optim_pdds_robust} being the parallel version.
}
\usage{
optim_pso        (objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, max_number_function_calls=NULL, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile = "ppso.pro", save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE)
optim_ppso       (objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, max_number_function_calls=NULL, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile = "ppso.pro", save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL,                                     nslaves = -1)
optim_ppso_robust(objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, max_number_function_calls=NULL, w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)), Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile = "ppso.pro", save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves = -1, working_dir_list=NULL)

optim_dds        (objective_function = sample_function, number_of_parameters = 2, number_of_particles =  1,                               max_number_function_calls= 500, r=0.2,                 abstol = -Inf, reltol = -Inf, max_wait_iterations=50,                                    parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)),                                                           lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile=  "dds.log",  projectfile = "dds.pro",  save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE)
optim_pdds_robust(objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40,                               max_number_function_calls= 500, r=0.2,                 abstol = -Inf, reltol = -Inf, max_wait_iterations=50,                                    parameter_bounds = cbind(rep(-1, number_of_parameters), rep(1, number_of_parameters)),                                                           lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile=  "dds.log",  projectfile = "dds.pro",  save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves = -1, working_dir_list=NULL,...)

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

  \item{max_number_function_calls}{
Abortion criteria: Maximum total number of calls to the objective function, ignored when set to \code{NULL}.
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

  \item{plot_progress}{
If set to TRUE, corresponds to repeated call of \code{\link{plot_optimization_process}} with default parameters in an interval of \code{save_interval}. If \code{plot_progress} is a named list, it is parsed as paraemters to \code{\link{plot_optimization_process}}.
}


  \item{tryCall}{
If set to \code{TRUE}, \code{objective_function} is executed using \code{try}, which yield error messages if \code{objective_function} fails on any slave. For very fast evaluating functions, this may increase evaluation time.
}



#  \item{...}{
#Further arguments passed to \code{objective_function}.
#}

for \code{optim_ppso*}:

  \item{w}{
Inertia constant, i.e. "weight" of the particles.
}
  \item{C1}{
Cognitive component. Weighting factor for the "personal" experience of each particle.
}
  \item{C2}{
Social component. Weighting factor for the "swarm" experience of each particle.
}

for \code{optim_pdds*}:

  \item{r}{
Neighbourhood size perturbation parameter . Default value: 0.2.
}


for the parallel versions \code{optim_ppso*}, \code{optim_pdds*}:
  \item{nslaves}{
number of rmpi slaves to spawn (default -1: as many as possible, requires package \code{Rmpi})
}
  \item{working_dir_list}{
 String matrix of working directories the slaves should change to. First column holds the hostname, the second column the directory.\\ 
 A hostname may be listed several times, if more than one slave is run on it. If a host is listed less times than it has slaves, the last entry will be recycled as necessary.\\ 
 The entry "default" denotes the directory all slaves of unlisted hosts will be change to. If missing, slaves on unlisted hosts will stay in the directory they are spawned in (don't ask me which this is).\\
 Beware: "~" apparently doesn't work, so rather use the full path. 
}

(didactical purpose, only for two-parameter search and fast objective function, non-parallelized. )
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
\item{optim_pso}{ serial computation for fast objective functions and single core machines}
\item{optim_ppso}{			 parallel computation: computationally demanding objective function and multi-core CPU or network slaves.       }
\item{optim_ppso_robust}{slightly slower than \code{optim_ppso}, but presumably more robust (to be tested), since it doesn't leave slaves in "listen" mode after crashes}
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
  Tolson, B. A., and C. A. Shoemaker (2007)
  Dynamically dimensioned search algorithm for computationally efficient watershed model calibration, \emph{Water Resour. Res.,}
  43, W01413, doi:10.1029/2005WR004723.  \url{http://www.agu.org/journals/wr/wr0701/2005WR004723/}
}
\author{
Till Francke <francke_at_uni-potsdam.de>
}
\note{
- still in development, use with care, all comments welcome -
}


\seealso{
 \code{\link{plot_optimization_process}}, \code{\link{sample_function}}
}

\examples{
 library(ppso)
#simple application (all file I/O disabled)
	result = optim_pso(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL)
	print (result)        #actual minimum -2 at (0,0) 
	result = optim_dds(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL)
	print (result)


#simple application with visualisation
	result = optim_pso(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL, do_plot="rgl")
	print (result)
	result = optim_dds(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL, do_plot="base")
	print (result)


#writing and resuming from project file
	projectfile	=tempfile()
	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=NULL, load_projectfile="no")		#start optimization, generate project file
	print (result)

	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=NULL, load_projectfile="yes")		#resume optimization from previous project file
	print (result)
	unlink(projectfile)	#delete project file

#visualisation of progress
	projectfile	=tempfile()
  logfile	=tempfile()
	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=logfile, load_projectfile="no", plot_progress=TRUE)		#start optimization, generate project file and log file
  unlink(c(projectfile,logfile))	#delete log file and project file


#parallel application
	result = optim_ppso(objective_function=rastrigin_function, nslaves=2, max_number_function_calls=200, projectfile=NULL, logfile=NULL)     #standard version
	print (result)

	result = optim_ppso_robust(objective_function=rastrigin_function, nslaves=2, max_number_function_calls=200, projectfile=NULL, logfile=NULL)    #"robust" version
	print (result)

	result = optim_pdds_robust(objective_function=rastrigin_function, nslaves=2, max_number_function_calls=200, projectfile=NULL, logfile=NULL)    #"robust" version
	print (result)
	
  working_dir_list=rbind( #specify working directories for four slaves
  c(host="default", wd="/tmp/defaultdir"),
  c(host="host1", wd="/home/me/firstthread"),
  c(host="host1", wd="/home/me/firstthread"),
  c(host="host2", wd="/home/me/onlyonethread")
  )
	result = optim_pdds_robust(objective_function=rastrigin_function, nslaves=4, working_dir_list=working_dir_list)    


}

\keyword{ ~optimization }
%\keyword{ ~parallel } 
%\keyword{ ~swarm } 
%\keyword{ ~particle } 
%\keyword{ ~dynamically } 
%\keyword{ ~dimensioned } 
%\keyword{ ~search } 


% __ONLY ONE__ keyword per line