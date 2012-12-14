\name{push_object}
\Rdversion{1.1}
\alias{push_object}
\title{
Pushing an object to the master.
}
\description{
push an object to the master.
}
\usage{
push_object (object_name, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object_name}{
name of object to be pushed to master.
}
\item{value}{
value of object to be set on the master.
}
}
\details{
  This function pushes an object to the master. This is not necessary for ordinary optimizations! However, it enables the slave to react specifically on the overall conditions of the optimization.
  Potentially useful objects are the variables listed and explained in globvars.R (see package source code).
 }
\value{
none.
}
\author{
Till Francke
}
\note{
  The pushed object is assigned in the local frame of the calling master function, or, if not found there, in the globvars environment.
}

\seealso{
 \code{\link{request_object}}
}

\examples{

\dontrun{
#objective function with additional master-slave interaction
interact_function = function(parms)
{
  library(ppso)
  X=push_object("X")   #push some object
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
