#internal function: empty message queues and close Rmpi session
close_mpi=function()                        
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{
  while(mpi.iprobe(mpi.any.source(),5))                #answer pending object requests from slaves (otherwise they cannot be shut down gracefully)
  {
      slave_message <- mpi.recv.Robj(mpi.any.source(), 5)
      slave_message_info <- mpi.get.sourcetag()
      slave_id <- slave_message_info[1]
     	mpi.send.Robj(obj="kill", dest=slave_id, tag=7)         #send kill signal
  }

  i=0
  if (mpi.comm.size()>0)        #somehow, the previous loops crashes the comm
  {
    while(mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #empty MPI queue if there is still something in there
    {
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        i=i+1 #count orphaned messages
    }
    if (i > sum(globvars$status==2)) warning(paste(i,"orphaned MPI-messages discarded."))
  }
  if (i < sum(globvars$status==2)) 
  {  
    warning(paste(sum(globvars$status==2)-i,"slave(s) may still be evaluating."))
    if (verbose_master) print(paste(sum(globvars$status==2)-i,"slave(s) may still be evaluating.")) 
  }
  
  if (globvars$mpi_mode=="loop")
  {  
    if (verbose_master){ print(paste(Sys.time()," ...sending slaves kill-tag...")); flush.console()}
    for (slave_id in 1:globvars$nslaves) mpi.send.Robj(obj="kill", dest=slave_id, tag=7)
  }
  #mpi.abort() #causes hangup
  
  if (!is.null(globvars$nslaves) &                                                    #there are slaves
        ((globvars$closed_slaves ==0) || (globvars$closed_slaves == globvars$nslaves))) #all (or none) are still available
#        & all(globvars$status==1)) 
  {
    if (verbose_master){ print(paste(Sys.time()," ...telling slaves to close log...")); flush.console()}
    if (verbose_slave) mpi.bcast.cmd(sink()) #terminate file logs of slaves - causes hangup for busy slaves

    if (verbose_master){ print(paste(Sys.time()," ...trying to close slaves...")); flush.console()}
    
    #mpi.close.Rslaves()    #causes hang, when slaves are still in mpi_loop (which they don't leave, because the kill-message has been deleted?)
    
    #mpi.quit()   #causes hang on orson
    if (verbose_master){ print(paste(Sys.time()," end of close_mpi.")); flush.console()}
  }  
  else
    if (interactive()) warning("Couldn't close some slaves. You may try by calling 'mpi.close.Rslaves()'. BEWARE: This may crash R, sorry for the inconvenience.")
}
