#internal function: prepare MPI connection and slaves

#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

#tag-codes for master-slave communications:
#0: task from master
#1: 
#2: results returning from slave
#3: slave sends good-bye message
#4: execution error from slave
#5: object request from slave or subsequent object transfer from master to slave
#6: object pushed from slave to master
#7: command slaves to abort everything

globvars$mpi_mode="loop" #"loop"    ("bcast" is probably obsolete, as it is less stable and slower)

prepare_mpi_cluster=function(nslaves, working_dir_list=NULL, verbose_slave=FALSE, ...)
{
  if (!is.loaded("mpi_initialize")) {         
  	if (!require("Rmpi")) stop("Package Rmpi not found. Install it or use serial version of this call (optim_pso or optim_dds).")
  }
 
  if (nslaves == -1) nslaves=mpi.universe.size() else  # Spawn as many slaves as possible
  if (nslaves > mpi.universe.size()) warning("Number of specified slaves exceeds number of available slaves.")

  if (mpi.comm.size()>0)
  {
    nslaves=mpi.comm.size()-1
	  print(paste(nslaves,"running slaves detected, no spawning."))
  } else {
 	  if (verbose_master) {print("spawning slaves..."); flush.console()}

    mpi.spawn.Rslaves(nslaves=nslaves)
  	print(paste(mpi.comm.size(),"slaves spawned."))
  }
  globvars$is_mpi = TRUE
  
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
  
   options(error=.Last)     #close rmpi on errors

  
  perform_task = function(params,slave_id) {
    
    # put all arguments into a list
    argl <- additional_args
    argl[[length(argl)+1]] <- params


      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": received message for slave",slave_id))

      if (!(mpi.comm.rank() %in% slave_id)) return(1) #only the adressed slaves should attend to the task
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": task received"))

      if (tryCall)
      {
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": calling objective function..."))
        #results=try(objective_function(params),silent=TRUE)  # call the objective function with the respective parameters, and create results (with error handling, slower)
        results=try(do.call(what = objective_function, args=argl),silent=TRUE)  # call the objective function with the respective parameters, and create results (with error handling, slower)
        
        if (verbose_slave)
    		{
    			print(paste(Sys.time(),"slave",mpi.comm.rank(),":  ...objective function evaluation completed"))
    			flush.console()
    		}

        if (is.numeric(results))
        {
          if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": returning results to master..."))
          mpi.send.Robj(results,0,2)      # Send the results back as a task_done slave_message            
        } else                             #an error occured during execution
        {
          if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": error during call of objective function, returning message to master. Error message:",as.character(results)))
          mpi.send.Robj(paste("(",Sys.info()["nodename"],"):",as.character(results)),0,4)    #return the error message, tagged as "error" (4)
        }        
      } else        #non-tryCall option
      {  
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": calling objective function..."))
        results=do.call(what = objective_function, args=argl)  # call the objective function with the respective parameters, and create results (without error handling, faster)
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...objective function evaluation completed"))   
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": returning results to master..."))
        mpi.send.Robj(results,0,2)      # Send the results back as a task_done slave_message            #ii isend doesn't work - why?
      }

      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...results returned, back to idle mode."))    
      return()
  } #end function "perform_task"
  
  additional_args <- list(...) #get list of arguments hidden in ...
  mpi.bcast.Robj2slave(additional_args)                   #send any additional arguments to slaves

  mpi.bcast.Robj2slave(verbose_slave)                   #send verbose-flags to slaves
  if (verbose_slave)
    mpi.bcast.cmd(sink(paste("slave",mpi.comm.rank(),sep="")))               #put output of slaves into files, if desired


  mpi.bcast.Robj2slave(tryCall)                     #tryCall-flag
  mpi.bcast.Robj2slave(objective_function)         #send objective function to slaves
  mpi.bcast.Robj2slave(perform_task)               #send activation function to slaves
  
  globvars$closed_slaves=0
  globvars$idle_slaves=1:nslaves

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
      globvars$closed_slaves=length(failed_dirchanges)         #count the slaves with failed dirchange as "closed"
      nslaves=nslaves-globvars$closed_slaves
      globvars$idle_slaves=globvars$idle_slaves[-failed_dirchanges]     #remove the slaves from list of ready slaves
    }
  }

  globvars$nslaves=nslaves
  
  if (globvars$mpi_mode=="loop")     #alternative slave-MPI-mode
  {
    if (verbose_master) print(paste(Sys.time(),"initiating slaves in loop-mode."))      
    mpi_receive_loop = function()
    {
      tag=0
      while (tag !=7 )
      {
          messge = mpi.recv.Robj(source=0, tag=mpi.any.tag())
          messge_info = mpi.get.sourcetag()
          tag      = messge_info[2]

          if ((tag == 0) )     #task message
            perform_task(params=messge, slave_id=mpi.comm.rank())        #do task
            
          if (((tag == 7) && (messge == "kill")) ||
              !is.null(globvars$kill_msg))     #kill-message
            mpi.send.Robj(obj="bye", dest=0, tag=3) #tag 3 demarks good-bye


      } 
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": leaving loop, good bye."))    
    }

    mpi.bcast.Robj2slave(mpi_receive_loop)               #transfer mpi-loop function to slaves
    mpi.bcast.cmd(mpi_receive_loop())                    #call     mpi-loop function on slaves
  }
  
  	
}


send_task = function(params, slave_id=slave_id)        #submit job to slave
#this function is used to transfer a task to a slave
#two different modi (globvars$mpi_mode="bcast" or "loop") for running the slaves are implemented, so there are two branches
{
  if (globvars$mpi_mode=="bcast")
   mpi.remote.exec(cmd=perform_task, params=params, slave_id=slave_id, ret=FALSE)        #submit job to slave by broadcasting to all
   else
   mpi.send.Robj(obj=params, dest=slave_id, tag=0) #send params to specified slave with tag=0 (this is a task)
}   


