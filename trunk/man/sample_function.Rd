\name{sample_function}
\Rdversion{1.1}
\alias{sample_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A sample function that can be optimized for 2 parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample_function(param_vec, fast = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param_vec}{
%%     ~~Describe \code{param_vec} here~~
}
  \item{fast}{
%%     ~~Describe \code{fast} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(param_vec,fast=FALSE)          #objective function to be minimized
{
  obj=0
  for (i in 1:NROW(param_vec))
    obj=obj+(cos(param_vec[i])+1)
  obj=-abs(obj)

 # r <- sqrt((param_vec[1]*4)^2+(param_vec[2]*4)^2) 
#  obj=-5 * sin(r)/r 

  if (!fast)
  {
    maxwaittime=0.1               #function to be run on slaves, with random delay
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
