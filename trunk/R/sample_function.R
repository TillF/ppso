sample_function <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
  obj=0
  for (i in 1:NROW(param_vec))
    obj=obj+(cos(param_vec[i])+1)
  obj=-abs(obj)


#  if (runif(1,0,1)>0.95)
#  {
#    cat("error in objective function simulated by slave ",mpi.comm.rank(), file=paste("slave",mpi.comm.rank(),".log",sep=""))
#    error("error in objective function simulated by slave ",mpi.comm.rank())
#  }

#  if (runif(1,0,1) > 0.8 && mpi.comm.rank()==1)
#    maxwaittime = 1          #simulate a very slow run on slave 1
    
  if (maxwaittime > 0)  #function to be run on slaves, with random delay
    Sys.sleep(maxwaittime)
  
  return(obj)
}

