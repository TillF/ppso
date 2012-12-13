\name{plot_optimization_progress}
\Rdversion{1.1}
\alias{plot_optimization_progress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualisation of optimization process
}
\description{
Visualisation of optimization process based on log files. \code{plot_optimization_progress} is also used by {\code{\link{optim_pso}}}  & Co. when used with \code{plot_optimization_progress=TRUE}.
}
\usage{
plot_optimization_progress (logfile="pso.log", projectfile="pso.pro", progress_plot_filename=NULL, goodness_plot_filename=NULL, cutoff_quantile=0.95, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logfile}{
log file as produced by \code{optim_pso} & Co.
}
  \item{projectfile}{
project file as produced by \code{optim_pso} & Co.
}
  \item{progress_plot_filename}{
optional file name to save the contents of the progress_plot window to. See details. 
}
  \item{goodness_plot_filename}{
optional file name to save the contents of the goodness_plot window to. See details. 
}
  \item{cutoff_quantile}{
upper quantile at which poorer goodness values are no more displayed (e.g. 0.95). This clarifies scatter plots.
}
  \item{verbose}{
print information on current optimum, the respective parameters and its function call number.
}
}

\details{
\code{plot_optimization_progress} can be used to inspect an ongoing optimization progress if its recording of log file and project file are enabled.
It produces two windows:
\itemize{
  \item{progress_plot}{
depicts the variation of each parameter versus the number of function evaluations. This helps checking which of the parameters are converging. Likewise, the development of the value of the objective function is illustrated. The best parameter set of each particle so far is marked with a red diamond.
}
  \item{goodness_plot}{
depicts the value of each parameter versus the objective function, giving clues to parameter sensitivity. Additionally, the distribution of the workload of the slaves and their time for function evaluations is displayed. Execution times less than a second are displayed as 0.5 second.
}
}
The plots are saved when \code{progress_plot_filename} or \code{goodness_plot_filename}, respectively, are specified. This is especially useful in conjunction with internal calls by {\code{\link{optim_pso}}} & Co. to supervise an optimization run. If the file names contain an extension supported by savePlot, this file format is used, otherwise a PostScript file is stored.
}

\value{
vector of optimum solution so far.
}

\author{
Till Francke
}

\seealso{\code{\link{optim_pso}}} 

\examples{
\dontrun{
	projectfile	=tempfile()
  logfile	=tempfile()
	result = optim_pso(objective_function=rastrigin_function, projectfile=projectfile, logfile=logfile, load_projectfile="no")		#start optimization, generate project file and log file
  plot_optimization_progress(logfile=logfile, projectfile=projectfile)
  unlink(c(projectfile,logfile))	#delete log file and project file
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~optimization }


