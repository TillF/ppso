\name{request_object}
\Rdversion{1.1}
\alias{request_object}
\title{
request_object
}
\description{
Request an object from the master.
}
\usage{
request_object (object_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object_name}{
name of object to be requested from master.
}
}
\details{
  This function requests an object from the master. This is not necessary for ordinary optimizations! However, it enables the slave to react specifically on the overall conditions of the optimization.\cr
  Potentially useful objects are the variables listed and explained in globvars.R (see package source code).
 }
\value{
If existing, returns the requested object, otherwise \code{NA}}

\author{
Till Francke
}
\note{
In serial mode (optim_*), this function basically looks up the search tree and the package environment for the requested objects. In parallel mode (optim_p*), mpi-communications are used.\cr
}

\seealso{
 \code{\link{push_object}}
}

\examples{

\dontrun{
#objective function with additional master-slave interaction
interact_function = function(parms)
{
  library(ppso)
  X=request_object("X")   #request some object
  print(paste("X:",X))
  flush.console()
  
  print(paste("pushing object..."))
  push_object(object_name="fitness_gbest", value=1000)     #push some object
  flush.console()
  
  print("returning results...")
  return(ackley_function(parms))
}

result = optim_pdds_robust(objective_function=interact_function, projectfile=NULL, logfile=NULL, load_projectfile="no",
                            max_number_function_calls=20, verbose=TRUE, tryCall=TRUE)
}
}							
							
							
							
\keyword{ ~optimization }
%\keyword{ ~example }% __ONLY ONE__ keyword per line
