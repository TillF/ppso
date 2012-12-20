\name{check_execution_timeout}
\Rdversion{1.1}
\alias{init_particles}
\alias{init_visualisation}
\alias{prepare_mpi_cluster}
\alias{update_tasklist_dds}
\alias{update_tasklist_pso}


\title{
Internal functions, not to be called by user directly
}
\description{
Internal functions,\\ not to be called by user directly
}
\usage{
check_execution_timeout(maxtries=10)
init_particles(lhc_init = FALSE)
init_visualisation()
prepare_mpi_cluster(nslaves, working_dir_list = NULL, verbose_slave = FALSE)
update_tasklist_dds(loop_counter = 1)
update_tasklist_pso()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{}
}

%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
None.}
%% \references{}

\author{
Till Francke
}
\note{
   This is an incomplete list of internal functions, which are not to be called by the package user.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~example }% __ONLY ONE__ keyword per line
