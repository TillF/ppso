sample_function2 <-
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

