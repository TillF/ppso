#test function after Ackley (1987), minimum value=-20-e
ackley_function <-
function(param_vec,maxwaittime=0)          #objective function to be minimized
{
  obj=0
  dimen=length(param_vec);

  sum1 = sum(param_vec^2)
  sum2 = sum(cos(2*pi*param_vec))

  obj = -20*exp(-0.2*(sum1/dimen)^0.5)-exp(sum2/dimen);


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


