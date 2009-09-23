#test function after Griewank (1981), minimum value=0
griewank_function <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
  dimen=length(param_vec);

  obj=sum(param_vec^2/4000)- prod(cos(param_vec/sqrt(1:dimen)))+1

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


