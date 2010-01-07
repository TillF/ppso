#internal function: check if any of the nodes exceeds the allowed execution time
check_execution_timeout_i = function()
{
  maxtries=3        #maximum number of attempts to run on a slave before it is no longer used
  
  nodes=unique(node_id[status==2])                         #check all nodes that are currently employed

  for (i in nodes)
  {
    current_particle = node_id==i & (status==2)

    current_time_execution=difftime(Sys.time(),computation_start[current_particle],units="secs")       #get the time this node has been processing the current call so far
 
    calls_so_far=execution_times$slave_id==i                                 #calls performed by this node so far
    if(sum(calls_so_far)==0)                                             #no completed calls so far
      if (nrow(execution_times)>=3*length(nodes))                           #have there been sufficiently many other calls?
        mean_execution_time=mean(execution_times$secs)                      #use overall mean as benchmark
      else      
        next                                                                #no calls so far, simply wait
    else
      mean_execution_time=mean(execution_times$secs[calls_so_far])    #calculate the mean execution time of this node so far

  
    if (current_time_execution > mean_execution_time*execution_timeout)     #current call takes too long?
    {
       node_id[current_particle] = 0          #reset particle
       status [current_particle] = 0  
       node_interruptions[i,"counter"] = node_interruptions[i,"counter"] +1   #increase counter of interruptions
       if (node_interruptions[i] > maxtries)
       {
        closed_slaves <- closed_slaves + 1
        warning(paste("Excluded node",i,"because of failing to produce results within",mean_execution_time,"s for",maxtries,"attempts."))
        node_interruptions[i,"status" ] = 2  #flag as "terminated permanently"
       } else
       node_interruptions[i,"status" ] = 1  #flag as "terminated once"
                              
    }
  }

  #"export" the variables
  assign("status",             status,           parent.frame())
  assign("node_id",            node_id,          parent.frame())
  assign("idle_slaves",       idle_slaves,       parent.frame())       
  assign("node_interruptions",node_interruptions,parent.frame())
  assign("closed_slaves",closed_slaves,          parent.frame())

  if (length(execution_times$secs)==0)
    return(0.1)
  else
    return(min(c(1,execution_times$secs)))
#return(0)     #rr
}
