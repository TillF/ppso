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
\code{optim_pso} minimizes a given function \code{objective_function} with regard to its parameters contained in the vector \code{par} towards a minimum value of f using Particle Swarm Optimization.
\code{optim_ppso_robust} is the parallelized versions (using multiple CPUs).
\code{optim_dds} minimizes using a Dynamically Dimensioned Search with  \code{optim_pdds_robust} being the parallel version.
}
\usage{
optim_pso        (objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, max_number_function_calls=500,
  w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters),
  rep(1, number_of_parameters)), initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile = "ppso.pro", save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, verbose=FALSE, ...)
optim_ppso_robust(objective_function = sample_function, number_of_parameters = 2, number_of_particles = 40, max_number_of_iterations = 5, max_number_function_calls=500,
  w = 1, C1 = 2, C2 = 2, abstol = -Inf, reltol = -Inf, max_wait_iterations = 50, wait_complete_iteration = FALSE, parameter_bounds = cbind(rep(-1, number_of_parameters),
  rep(1, number_of_parameters)), initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = "ppso.log", projectfile = "ppso.pro", save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves = -1, working_dir_list=NULL, execution_timeout=NULL, maxtries=10, verbose=FALSE, ...)

optim_dds        (objective_function = sample_function, number_of_parameters = 2, number_of_particles =  1,                               max_number_function_calls= 500,
  r=0.2,                 abstol = -Inf, reltol = -Inf, max_wait_iterations=50,                                    parameter_bounds = cbind(rep(-1, number_of_parameters),
  rep(1, number_of_parameters)), initial_estimates=NULL, part_xchange=2,                                                          lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile=  "dds.log",  projectfile = "dds.pro",  save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, verbose=FALSE, ...)
optim_pdds_robust(objective_function = sample_function, number_of_parameters = 2, number_of_particles = 1,                               max_number_function_calls= 500,
  r=0.2,                 abstol = -Inf, reltol = -Inf, max_wait_iterations=50,                                    parameter_bounds = cbind(rep(-1, number_of_parameters),
  rep(1, number_of_parameters)), initial_estimates=NULL, part_xchange=2,                                                          lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile=  "dds.log",  projectfile = "dds.pro",  save_interval = ceiling(number_of_particles/4), load_projectfile = "try", break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves = -1, working_dir_list=NULL, execution_timeout=NULL, maxtries=10, verbose=FALSE, ...)

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
Number of particles used in the Particle Swarm Optimization / number of parallel threads in DDS.
}
  \item{max_number_of_iterations}{
Abortion criterion: Maximum number of iterations (i.e. update of particle positions) that are used.
}

  \item{max_number_function_calls}{
Abortion criterion: Maximum total number of calls to the objective function, ignored when set to \code{NULL}. In non-parallel mode, slightly more function calls may actually be made to finish the last iteration.
when resuming from a previous run, the previously made function calls also count, i.e. only the pending number of calls to reach \code{max_number_function_calls} will be made. In contrast, if \code{max_number_function_calls} is negative, -\code{max_number_function_calls} MORE calls will be made, additionally to those of the previous run. 
}

  \item{abstol}{
Abortion criterion: minimum absolute improvement between iterations  (default: -Inf)
}
  \item{reltol}{
Abortion criterion: minimum absolute relative improvement between iterations (default: -Inf)
}
  \item{max_wait_iterations}{
Number of iterations, within these an improvement of the above described (\code{abstol, reltol}) quality has to be achieved before abortion
}
  \item{wait_complete_iteration}{
TRUE:  Update particle velocities after function \code{objective_function} has been evaluated for ALL particles.
FALSE: Update particle velocities after each completed evaluation of function \code{objective_function} (suggested in parallel mode).
}
  \item{parameter_bounds}{
(row-named) matrix containing lower (first column) and upper (second column) boundary for each parameter. See details.
}
  \item{initial_estimates}{
(row-named) matrix containing columns of initial estimates (column-wise). Pending computations from a resumed project file are used first, then all \code{initial_estimates} are used.
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
Name of a file file that will cause the termination of the optimization when it is encountered. Useful for interrupting a running optimization gracefully. At start-up, any existing file of this name will be deleted.
}

  \item{plot_progress}{
If set to TRUE, corresponds to repeated call of \code{\link{plot_optimization_progress}} with default parameters in an interval of \code{save_interval}. If \code{plot_progress} is a named list, it is parsed as paraemters to \code{\link{plot_optimization_progress}}.
}

  \item{tryCall}{
If set to \code{TRUE}, \code{objective_function} is executed using \code{try}, which yield error messages if \code{objective_function} fails on any slave. The optimization retries, but terminates when an error is produced a second time for the same parameter set. For very fast evaluating functions, this may increase evaluation time.
} 

\item{...}{
additional arguments passed to \code{objective_function}
}


for \code{optim_*pso*}:

  \item{w}{
Inertia constant, i.e. "weight" of the particles.
}
  \item{C1}{
Cognitive component. Weighting factor for the "personal" experience of each particle.
}
  \item{C2}{
Social component. Weighting factor for the "swarm" experience of each particle.
}

for \code{optim_*dds*}:

  \item{r}{
Neighbourhood size perturbation parameter . Default value: 0.2.
}
\item{part_xchange}{
Relevant for \code{number_of_particles} > 1: mode how DDS particles (i.e. parallel DDS-threads) communicate with each other between iterations:
\itemize{
\item{0:} {no communication/relocation between particles}
\item{1:} {relocate/update particle particle that is worst in both objective function AND futile iterations}
\item{2:} {relocate all but the best particle}
\item{3:} {relocate particle that is worst in futile iterations (but not the global best) and set to best improving particle}
}
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
		\item{execution_timeout}{
		If set to a factor greater than 1, a slave whose execution time has exceeded (\code{execution_timeout} times the mean of its prior calls) for more than three times is no longer used. This also applies when the slave produces its first result very slow when compared to the other slaves. 
		}
		\item{maxtries}{
		Number of times a slave may exceed \code{execution_timeout} in a row before being excluded from further tasks. Ignored, when \code{execution_timeout} is not set.
		}
		\item{verbose}{
		   generate screen output documenting the message passing}
		   \itemize{
			   \item{"master"}{master generates output}
			   \item{"slaves"}{ slaves generate output (into text files named "slave<nn>")}
			   \item{TRUE}{ master and slaves generate output}
			   }


for didactical purpose (only for two-parameter search and fast objective function, non-parallelized):
  \item{do_plot}{
enable 3D-plot of response surface and search progress.
"base": use basic 3D-plotting functions.
"rgl":  use moveable rgl-plotting (requires package \code{rgl})
}
  \item{wait_for_keystroke}{
waiting for keystroke between iterations (e.g. for plotting). Option for changing into debug mode.
}
}


\details{
If \code{parameter_bounds} (dominant) or \code{initial_estimates} has named rows, these names are used in vector passed to the the call to \code{objective_function}.

recommended use:
\itemize{
\item{optim_pso, optim_dds}{ serial computation for fast objective functions and single core machines}
\item{optim_ppso_robust, optim_pdds_robust}{			 parallel computation: computationally demanding objective function and multi-core CPU or network slaves.       }
}
concerning reproducibility: setting the seed using \code{set.seed()} should make results reproducable. However, in the parallel versions, the randomness introduced by the order of the results returned from the slaves cannot be controlled, so this may limit reproducibility.

}

\value{
The functions return a list with the elements
  \item{par }{parameters (ie location) of minimum found}
  \item{value }{value of objective function at minimum}
  \item{function_calls }{total number of function calls performed}
  \item{break_flag }{criterium that caused the termination of the algorithm (\code{c("max number of function calls reached", "max iterations reached", "converged", "all slaves gone", "user interrupt")})}
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

Parallelization is a great thing, especially when it works. On its dark side, you just parallelize errors and multiply debugging effort. Enabling \code{verbose=TRUE} and  \code{tryCall=TRUE} may help to find the error.\cr
Master or slave crashes can cause the entire session to stall, so restarting R may be necessary. Under Windows, orphaned mpi-sessions and slaves may have to be killed manually to prevent stalling the entire system: Use "Process Explorer" (free) and kill the smpd.exe-tree and mpiexec.exe.
}


\seealso{
 \code{\link{plot_optimization_progress}}, \code{\link{sample_function}}
}

\examples{
 library(ppso)
#simple application (all file I/O disabled)
	result = optim_pso(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL)
	print (result)        #actual minimum -2 at (0,0) 
	result = optim_dds(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL)
	print (result)


\dontrun{
#simple application with visualisation
	result = optim_pso(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL, do_plot="rgl")
	print (result)
	result = optim_dds(objective_function=rastrigin_function, projectfile=NULL, logfile=NULL, do_plot="base")
	print (result)
}

#writing and resuming from project file
	projectfile	=tempfile()
	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=NULL, load_projectfile="no")		#start optimization, generate project file
	print (result)

	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=NULL, load_projectfile="yes")		#resume optimization from previous project file
	print (result)
	unlink(projectfile)	#delete project file

\dontrun{
#visualisation of progress
	projectfile	=tempfile()
  logfile	=tempfile()
	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=logfile, load_projectfile="no", plot_progress=TRUE)		#start optimization, generate project file and log file
  unlink(c(projectfile,logfile))	#delete log file and project file



#parallel application
	result = optim_ppso_robust(objective_function=rastrigin_function, nslaves=2, max_number_function_calls=200, projectfile=NULL, logfile=NULL)
	print (result)

	result = optim_pdds_robust(objective_function=rastrigin_function, nslaves=2, max_number_function_calls=200, projectfile=NULL, logfile=NULL)
	print (result)
	
  working_dir_list=rbind( #specify working directories for four slaves
  c(host="default", wd="/tmp/defaultdir"),
  c(host="host1", wd="/home/me/firstthread"),
  c(host="host1", wd="/home/me/secondthread"),
  c(host="host2", wd="/home/me2/onlyonethread")
  )
	result = optim_pdds_robust(objective_function=rastrigin_function, nslaves=4, working_dir_list=working_dir_list)    
}

}

\keyword{ ~optimization }
%\keyword{ ~parallel } 
%\keyword{ ~swarm } 
%\keyword{ ~particle } 
%\keyword{ ~dynamically } 
%\keyword{ ~dimensioned } 
%\keyword{ ~search } 
