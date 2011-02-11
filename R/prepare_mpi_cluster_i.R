#internal function: prepare MPI connection and slaves

#note: this function reads and writes to non-local variables (i.e. varaibles declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_p*  (clumsy, but I don't know better)

prepare_mpi_cluster_i=function(nslaves=nslaves, working_dir_list=NULL, verbose=FALSE)
{
  if (!is.loaded("mpi_initialize")) {         
  library("Rmpi")
  }
 
  if (nslaves==-1) nslaves=mpi.universe.size() else  # Spawn as many slaves as possible
  if (nslaves>mpi.universe.size()) warning("Number of specified slaves exceeds number of available slaves.")

  mpi.spawn.Rslaves(nslaves=nslaves)
  
  while(mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #empty MPI queue if there is still something in there
    slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())


     .Last <- function(){
      if (is.loaded("mpi_initialize")){
          if (mpi.comm.size(1) > 0){
              #print("Please use mpi.close.Rslaves() to close slaves.")
              mpi.close.Rslaves()
          }
          #print("Please use mpi.quit() to quit R")
          #.Call("mpi_finalize")
      }
    }
  
   print(paste(mpi.comm.size()-1,"slaves spawned"))
   options(error=.Last)     #close rmpi on errors

  
  #options(error=.Last)     #close rmpi on errors

  perform_task <- function(params,slave_id,tryCall=FALSE) {

      if (verbose) print(paste(Sys.time(),"received message for slave",slave_id))
      #write.table(file=paste("slave",slave_id),"contacted")

      if (!(mpi.comm.rank() %in% slave_id)) return(1) #only the adressed slaves should attend the task
      if (verbose) print(paste(Sys.time(),"task received"))


      #write.table(file=paste("slave",slave_id),"task received")
      if (tryCall)
      {
        if (verbose) print(paste(Sys.time(),"calling objective function..."))
        results=try(objective_function(params),silent=TRUE)  # call the objective function with the respective parameters, and create results (with error handling, slower)
        if (verbose) print(paste(Sys.time(),"...objective function evaluation completed"))

        if (!is.numeric(results))                      #an error occured during execution
        {
          if (verbose) print(paste(Sys.time(),"returning results to master..."))
          mpi.send.Robj(paste("(",Sys.info()["nodename"],"):",as.character(results)),0,4)    #return the error message, tagged as "error" (4)
          if (verbose) print(paste(Sys.time(),"...results returned, back to idle mode."))
          return()
        }        
      }
      else
      {  
        if (verbose) print(paste(Sys.time(),"calling objective function..."))
        results=objective_function(params)  # call the objective function with the respective parameters, and create results (without error handling, faster)
        if (verbose) print(paste(Sys.time(),"...objective function evaluation completed"))   
        if (verbose) print(paste(Sys.time(),"returning results to master..."))
        mpi.send.Robj(results,0,2)      # Send the results back as a task_done slave_message            #ii isend doesn't work - why?
        if (verbose) print(paste(Sys.time(),"...results returned, back to idle mode."))    
        return()
      }
        
  }
  

  mpi.bcast.Robj2slave(objective_function)         #send objective function to slaves
  mpi.bcast.Robj2slave(verbose)                   #send verbose-flags to slaves
  mpi.bcast.Robj2slave(perform_task)               #send activation function to slaves

  closed_slaves=0
  nslaves=nslaves
  idle_slaves=1:nslaves

  if (!is.null(working_dir_list))    #change working directories of slaves, if specified
  {
    slave_hosts=capture.output(slave.hostinfo())
    slave_hosts=sub(" ","",sub("[^:]*: ","",slave_hosts[-1]))  #extract hostnames of slaves, discard master

    wd_to_be_set=array("",length(slave_hosts))                                       #this will be a list with a default directory to each slave
    
    default_dir=working_dir_list[which(working_dir_list[,"host"]=="default"),"wd"]  #get name of default_dir
    if (length(default_dir)==0) default_dir="./"                                    #if default dir not specified, use "./"
     
    for (i in 1:length(slave_hosts))
    {
      entry=which(working_dir_list[,"host"]==slave_hosts[i])
      if (length(entry)==0)
        wd_to_be_set[i]=default_dir                                                 #set default dir
      else
        wd_to_be_set[i]=working_dir_list[entry[1],"wd"]                             #set specified dir
      if (length(entry)>1) working_dir_list=working_dir_list[-entry[1],]            #the current host was listed more than once in the list, remove the first entry  
    }    
    
    mpi.bcast.Robj2slave(wd_to_be_set)            #sent the directory list to all slaves
    res=mpi.remote.exec(setwd(wd_to_be_set[mpi.comm.rank()]))
    #failed_dirchanges= (as.character(mpi.remote.exec(getwd()))!=wd_to_be_set)     #check which of the dirchanges did not succeed
    
    failed_dirchanges=NULL
    for (i in 1:length(res))
      if (!is.null(attr(res[[i]],"class")))
        failed_dirchanges=c(failed_dirchanges,i)
    
    if (any(failed_dirchanges))
    {
      warning(paste("The following slave(s) could not change into the respective directory and remain unused:\n",
        paste(t(cbind(slave_hosts," : ",wd_to_be_set,"\n ")[failed_dirchanges,]),collapse="")))
      closed_slaves=length(failed_dirchanges)         #count the slaves woth failed dirchange as "closed"
      nslaves=nslaves-closed_slaves
      idle_slaves=idle_slaves[-failed_dirchanges]     #remove the slaves from list of ready slaves
    }
  }

 
  assign("closed_slaves",closed_slaves,parent.frame())        #set globals
  nslaves = mpi.comm.size()-1
  assign("nslaves",         nslaves,   parent.frame())        #  
  assign("idle_slaves", idle_slaves,   parent.frame())        #
}


