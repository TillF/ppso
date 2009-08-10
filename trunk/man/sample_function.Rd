\name{sample_function}
\Rdversion{1.1}
\alias{sample_function}
\alias{sample_function2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample function that can be optimized for 2 parameters
}
\description{
Just demo functions for testin purposes.
}
\usage{
sample_function(param_vec, fast = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param_vec}{
parameter vector
}
  \item{fast}{
when set TRUE, the function is evaluated at maximum speed, otherwise a random delay of up to maxwaittime (=0.1 s) is introduced. Useful for simulating heterogenous clusters.   
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
Till Francke
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

##sample_function is currently defined as
function(param_vec,fast=FALSE)          #objective function to be minimized
{
  obj=0
  for (i in 1:NROW(param_vec))
    obj=obj+(cos(param_vec[i])+1)
  obj=-abs(obj)
 
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



## sample_function2 is currently defined as
function(param_vec,fast=FALSE)          #objective function to be minimized
{
  r <- sqrt((param_vec[1]*4)^2+(param_vec[2]*4)^2) 
  obj=-5 * sin(r)/r 

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
\keyword{ ~optimization }
\keyword{ ~example }% __ONLY ONE__ keyword per line
