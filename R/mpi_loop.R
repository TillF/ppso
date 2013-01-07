mpi_loop = function(init_search)   #loop in which the master coordinates slave action
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{

  while ((globvars$closed_slaves < globvars$nslaves) & is.null(globvars$break_flag))
  {
        if (verbose_master) flush.console()
        
        if (init_search)
        {if (all(globvars$status==1)) break} else     #for pre-search, no update of tasklist is needed but exit when all tasks are done
          update_tasklist_dds()    #update particle positions based on available results

        tobecomputed=globvars$status==0
        
        sort_index=sort(globvars$slave_status[globvars$idle_slaves, "counter"], index.return = TRUE)$ix  #sort idle slaves by number of interruptions
        globvars$idle_slaves = globvars$idle_slaves[sort_index]
  
        while ((length(globvars$idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
        {
           if (init_search)
           current_particle=which(tobecomputed)[1] else
           {
              current_particle=which.min(globvars$function_calls[tobecomputed])   #treat particles with low number of iterations first
              current_particle=which(tobecomputed)[current_particle[1]]     #choose the first entry
           } 
             
          slave_id=globvars$idle_slaves[1]                     #get free slave        
          if (verbose_master) print(paste(Sys.time()," ...sending task to slave",slave_id))  
          mpi.remote.exec(cmd=perform_task,params=globvars$X[current_particle,],tryCall=tryCall,slave_id=slave_id,ret=FALSE)        #submit job to slave
          if (verbose_master) print(paste(Sys.time()," ...task sent"))  
          globvars$idle_slaves=globvars$idle_slaves[-1]                         #remove this slave from list
          globvars$status            [current_particle]=2               #mark this particle as "in progress"
          globvars$node_id           [current_particle]=slave_id        #store slave_id of this task
          globvars$computation_start [current_particle]=Sys.time()      #store time of start of this computation

          if (!init_search)
          {
            if (verbose_master) print(paste(Sys.time()," ...updating task list"))  
            update_tasklist_dds()   #update particle speeds and positions based on available results
          }
          tobecomputed=globvars$status==0
        } 
  
          if (!is.null(globvars$break_flag) | all(globvars$status==1)) break
  
          sleeptime=0.1
          output_time=Sys.time()
          if (verbose_master) print(paste(Sys.time()," ...wait for messages from slaves..."))          
          while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()) &&              #wait till there is a MPI-message or...
                      !(any(globvars$status==0) && (length(globvars$idle_slaves)>0) ) )   #...there is something to do and free slaves are available             
          {
              if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
              {
                if (verbose_master) print("...break_file detected, aborting...")
                globvars$break_flag="user interrupt"   
                return()
              }  
              if (!is.null(globvars$execution_timeout)) sleeptime=check_execution_timeout()
              Sys.sleep(sleeptime) 
              if (verbose_master)
              {
                 if (difftime(Sys.time(), output_time, units="sec") > 10)
                 {
                  print(" ...still waiting...")
                  output_time=Sys.time()
                 }
              }
          }        
        if (!mpi.iprobe(mpi.any.source(),mpi.any.tag())) next #the previous loop was broken because of timeout, re-issue task
        if (verbose_master) {print(paste(Sys.time()," ...message detected")); flush.console()}
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        slave_message_info <- mpi.get.sourcetag()
        slave_id <- slave_message_info[1]
        tag      <- slave_message_info[2]
        if (verbose_master) print(paste(Sys.time()," ...message received from slave",slave_id, "(tag",tag,")"))          
#         browser()
        
        if (tag == 2) {      #retrieve result
          if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"sent results"))  
          current_particle =which(globvars$node_id==slave_id & globvars$status==2)           #find which particle this result belongs to
          if (length(current_particle) > 1)                                #rr   shouldn't occur
            print(paste("strange, slave",slave_id,"returned an ambiguous result (main search):",slave_message))
           else    
          if (length(current_particle) ==0)                                #
          {
            if (verbose_master) 
				if (globvars$slave_status[slave_id,"timeouts_in_row"] == 0)        #rr shouldn't occur
				   print(paste("strange, slave",slave_id,"returned an unrequested result (main search):",slave_message)) else
				   print(paste("slave",slave_id,"returned overdue result"))
          } else
          {
            globvars$fitness_X [current_particle] = slave_message
            globvars$status    [current_particle] =1      #mark as "finished"
            globvars$slave_status[slave_id,"timeouts_in_row"] = 0             #this result was in time, so reset the counter of consecutive timeouts 
            if (!is.null(globvars$execution_timeout)) globvars$execution_times = rbind(globvars$execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),globvars$computation_start[current_particle],units="sec"))))   #monitor execution times

            if (init_search)
            {
              if (!is.null(globvars$execution_timeout)) globvars$execution_times = rbind(globvars$execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),globvars$computation_start[current_particle],units="sec"))))   #monitor execution times
              if (!is.null(logfile))  write.table(file = logfile, cbind(format(globvars$computation_start[current_particle],"%Y-%m-%d %H:%M:%S"), matrix(globvars$X[current_particle,],ncol=number_of_parameters)  , globvars$fitness_X[current_particle], 
              globvars$node_id[current_particle]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)          
            }
            else
            globvars$function_calls[current_particle] =globvars$function_calls[current_particle]+1        #increase iteration counter
          }
          if (globvars$slave_status[slave_id,"timeouts_in_row" ] < maxtries) #if this slave failed too often in a row, don't use it any longer
            globvars$idle_slaves=c(globvars$idle_slaves,slave_id)
        } else if (tag == 3) {    # A slave has closed down.
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has gone away"))  
            globvars$closed_slaves = globvars$closed_slaves + 1
        } else if (tag == 4) {    # error occured during the execution of the objective function
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has produced an error:",slave_message))
  		    globvars$break_flag=paste("Abort, slave",slave_id,":",slave_message)
          globvars$closed_slaves=globvars$nslaves			
        } else if (tag == 5) {    #object request from slave 
    			slave_message_info <- mpi.get.sourcetag() 
    			 if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"requested object",slave_message))
          if(exists(x=slave_message))
    			  obj=get(slave_message,pos=parent.frame(), inherits=FALSE) else 
  			  if(exists(x=slave_message, where=globvars))
    			  obj=get(slave_message,pos=globvars) else 
            #return requested object, if existing, otherwise NA
    			  obj=NA
    			mpi.send.Robj(obj=obj, dest=slave_id, tag=5) #return object to slave
    			 if (verbose_master) print(paste(Sys.time()," ...object",slave_message, "sent to slave",slave_id))
		    } else if (tag == 6) {    #object pushed from slave 
    			slave_message_info <- mpi.get.sourcetag() 
    			object_name = attr(slave_message, "object_name") #name of object
   			  if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"pushed object",object_name,"=",slave_message))
#    			browser()
          if (!any(globvars$node_id==slave_id & globvars$status==2))
          {
             if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"is overdue, push ignored"))
             next
          }
          for (env in search())
            if(exists(x=object_name, where=env))                 #try to set this object in an environment, where it already exists
            {
              assign(x=object_name, value=slave_message, pos=env)
              if (verbose_master) print(paste(Sys.time()," ...object",object_name, "found in and assigned to",env))
              break
            }  
          if (!exists(x=object_name, inherits=TRUE))   #check if the object has not been found somewhere
          {
              assign(x=object_name, value=slave_message, pos=globvars)          
              if (verbose_master) print(paste(Sys.time()," ...object",object_name, "not found, assigned to globvars"))   			     
          }
	      } else  {    #unknown tag
   			 if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"sent unknown tag ",tag,", ignored"))
		    }  
  } 

  return()
}