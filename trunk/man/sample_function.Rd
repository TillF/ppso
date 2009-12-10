\name{sample_function}
\Rdversion{1.1}
\alias{sample_function}
\alias{sample_function2}
\alias{ackley_function}
\alias{griewank_function}
\alias{rastrigin_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample functions for testing purposes
}
\description{
Just demo functions that can be optimized for testing purposes.
}
\usage{
sample_function   (param_vec, maxwaittime = 0)
sample_function2  (param_vec, maxwaittime = 0)
ackley_function   (param_vec, maxwaittime = 0)
griewank_function (param_vec, maxwaittime = 0)
rastrigin_function(param_vec, maxwaittime = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param_vec}{
parameter vector
}
  \item{maxwaittime}{
when not set, the function is evaluated at maximum speed, otherwise a random delay of up to maxwaittime is introduced. Useful for simulating heterogenous clusters.   
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
Function value as defined by the respective function.}
\references{
  Ackley, D. H. (1987)
  A Connectionist Machine for Genetic Hillclimbing, Kluwer, Norwell, MA.
  
  Griewank, A. O. (1981)
  Generalized descent for global optimization,  \emph{Journal of Optimization Theory and Applications}, 34, 11-39.
  
  Rastrigin, L. A. (1974)
  Systems of Extremal Control, Nauka, Moscow.  \emph{in Russian}.
}

\author{
Till Francke
}
\note{
Minimum values:
sample_function:  -2*dimensions
sample_function2: -5
Ackley:	-20-e
Griewank: 0
Rastrigin: -dimensions
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }

\examples{

##sample_function is currently defined as
sample_function <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
  obj=0
  for (i in 1:NROW(param_vec))
    obj=obj+(cos(param_vec[i])+1)
  obj=-abs(obj)

 # r <- sqrt((param_vec[1]*4)^2+(param_vec[2]*4)^2) 
#  obj=-5 * sin(r)/r 

#  if (runif(1,0,1)>0.95)
#  {
#    cat("error in objective function simulated by slave ",mpi.comm.rank(), file=paste("slave",mpi.comm.rank(),".log",sep=""))
#    error("error in objective function simulated by slave ",mpi.comm.rank())
#  }

  if (maxwaittime>0)  #function to be run on slaves, with random delay
  {
    starttime=Sys.time()
    waittime= runif(1,0,maxwaittime)
    while (as.numeric(Sys.time()-starttime)<waittime)
    {
    }
  }
  
  return(obj)
}



## sample_function2 is currently defined as
sample_function2 <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
  r <- max(1e-320,sqrt((param_vec[1]*4)^2+(param_vec[2]*4)^2)) 
  obj=-5 * sin(r)/r 

  if (maxwaittime>0)  #function to be run on slaves, with random delay
  {
    starttime=Sys.time()
    waittime= runif(1,0,maxwaittime)
    while (as.numeric(Sys.time()-starttime)<waittime)
    {
    }
  }
  return(obj)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~optimization }
%\keyword{ ~example }% __ONLY ONE__ keyword per line
