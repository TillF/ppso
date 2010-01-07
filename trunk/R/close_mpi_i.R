#internal function: empty message queues and close Rmpi session
close_mpi_i=function()                        
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_* to allow accessing the same globals  (clumsy, but I don't know better)

{
  i=0
  while(mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #empty MPI queue if there is still something in there
  {
      slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
      i=i+1
  }
  if (i > sum(status==2)) warning(paste(i,"orphaned MPI-messages discarded."))
  if (i < sum(status==2)) warning(paste(sum(status==2)-i,"slave(s) may still be evaluating."))
  mpi.bcast.cmd(mpi.exit)
  if (!is.null(nslaves)) mpi.close.Rslaves()    #close cluster, if enabled
}
