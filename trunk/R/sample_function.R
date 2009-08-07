sample_function <-
function(param_vec,fast=TRUE)          #objective function to be minimized
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

