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

