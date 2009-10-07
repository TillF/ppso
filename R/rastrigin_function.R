#test function after Rastrigin (1974), minimum value=-dimen
rastrigin_function <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
#  dimen=length(param_vec);

  obj=sum(param_vec^2-cos(2*pi*param_vec))

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


