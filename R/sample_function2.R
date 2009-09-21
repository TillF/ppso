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

