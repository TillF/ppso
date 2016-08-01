#internal function: check if any of the nodes exceeds the allowed execution time
check_execution_timeout = function(verbose_master)
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
{
  #maxtries=10        #maximum number of attempts to run on a slave before it is no longer used
  
  #verbose_master = globvars$verbose_master #get global value

  nodes=unique(globvars$node_id[globvars$status==2])                         #check all nodes that are currently employed
  min_exp_runtime = Inf #minimum expected runtime of all active slaves
  
  for (i in nodes)
  {
    
    current_particle = globvars$node_id==i & (globvars$status==2)

    current_time_execution=difftime(Sys.time(),globvars$computation_start[current_particle],units="secs")       #get the time this node has been processing the current call so far
    
    calls_so_far = globvars$execution_times$slave_id==i                                 #calls performed by this node so far
    if(sum(calls_so_far)==0)                                             #no completed calls so far
      if (nrow(globvars$execution_times)>=3*length(nodes))                           #have there been sufficiently many other calls?
        mean_execution_time=mean(globvars$execution_times$secs)                      #use overall mean as benchmark
      else      
        next                                                                #no calls so far, simply wait
    else
      mean_execution_time = mean(globvars$execution_times$secs[calls_so_far])    #calculate the mean execution time of this node so far

    exp_runtime = max(0, mean_execution_time - current_time_execution) #expected runtime for this slave
    min_exp_runtime = min (exp_runtime, min_exp_runtime)       #minimum of all expected runtimes
    
  
    if (current_time_execution > mean_execution_time*globvars$execution_timeout)     #current call takes too long?
    {
       if (verbose_master) 
       {print(paste0(Sys.time(),": slave ",i," current execution time (", current_time_execution,") greater than", globvars$execution_timeout, "* benchmark (",mean_execution_time,"), resetting...")); flush.console()}

       globvars$node_id[current_particle] = 0          #reset particle
       globvars$status [current_particle] = 0 
       globvars$recent_error [current_particle] = i    #mark this particle as having produced an error recently
       
       globvars$slave_status[i,"counter"] = globvars$slave_status[i,"counter"] +1   #increase counter of total interruptions
       globvars$slave_status[i,"timeouts_in_row"] = globvars$slave_status[i,"timeouts_in_row"] +1   #increase counter of consecutive interruptions
       
       if (globvars$slave_status[i,"timeouts_in_row"] > maxtries)
       {
        globvars$closed_slaves = globvars$closed_slaves + 1
        warning(paste("Permanently excluded slave",i,"because of failing to produce results within",mean_execution_time*globvars$execution_timeout,"s for",maxtries,"attempts."))
       }                               
    } else
    if (verbose_master) {print(paste0(Sys.time()," ...current execution time (", current_time_execution,") less than ", globvars$execution_timeout, "* benchmark (",mean_execution_time,")")); flush.console()}
    
  }


  if (length(globvars$execution_times$secs)==0)
    return(0.1) #no measured execution times so far
  else
    return(min(30, min_exp_runtime)) #wait the minimum expected time until a slave finishes, but not more than 30 s
#return(0)     #rr
}
