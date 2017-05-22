mpi_loop = function(init_search=FALSE, method="dds", ...)   #loop in which the master coordinates slave action
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{
  verbose_master = globvars$verbose_master #get global value
  
if (method=="dds") update_tasklist= update_tasklist_dds else
                   update_tasklist= update_tasklist_pso 


  while ((globvars$closed_slaves < globvars$nslaves) & is.null(globvars$break_flag))
{
        if (verbose_master) flush.console()
     
       if (init_search)
        {if (all(globvars$status==1)) break} else     #for pre-search, no update of tasklist is needed but exit when all tasks are done
          update_tasklist()    #update particle positions based on available results

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
          if (verbose_master) {print(paste(Sys.time()," ...sending task to slave",slave_id)); flush.console()}  
            send_task(params=globvars$X[current_particle,], slave_id=slave_id)        #send job to slave
          if (verbose_master) {print(paste(Sys.time()," ...task sent")); flush.console()}  
            globvars$idle_slaves=globvars$idle_slaves[-1]                         #remove this slave from list
            globvars$status            [current_particle]=2               #mark this particle as "in progress"
            globvars$node_id           [current_particle]=slave_id        #store slave_id of this task
            globvars$slave_status[slave_id,"task_particle"]         =current_particle #store current task belonging to this slave
            globvars$computation_start [current_particle]=Sys.time()      #store time of start of this computation
             
         if (!init_search)
          {
            if (verbose_master) {print(paste(Sys.time()," ...updating task list")); flush.console()}  
            update_tasklist()   #update particle speeds and positions based on available results
          }
          tobecomputed=globvars$status==0
      } 
          if (!is.null(globvars$break_flag) | all(globvars$status==1)) break
  
          sleeptime=0.1
          output_time=Sys.time()
      
          if (wait_for_keystroke && (!exists("ch_mpi", where=globvars) || globvars$ch_mpi!="c")) 
          {
            if (is.numeric(globvars$ch_mpi) && (globvars$ch_mpi > 0)) globvars$ch_mpi=globvars$ch_mpi-1 else
            {
              print("MPI: ENTER to proceed, 'b'+ENTER for debug mode, <n> to skip <n> times, 'c'+ENTER to continue till end")
              globvars$ch_mpi=readline() 
            }  
            if (globvars$ch_mpi=="b")	browser() else
            if (globvars$ch_mpi!="c") globvars$ch_mpi = strtoi(globvars$ch_mpi) #counter
            if (is.na(globvars$ch_mpi)) globvars$ch_mpi=""
          }  
            

          if (verbose_master) {print(paste(Sys.time()," ...wait for messages from slaves...")); flush.console()}          
#          browser()
          while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()) &&              #wait till there is a MPI-message or...
                      !(any(globvars$status==0) && (length(globvars$idle_slaves)>0) ) )   #...there is something to do and free slaves are available             
          {
              if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
              {
                if (verbose_master) {print("...break_file detected, aborting..."); flush.console()}
                globvars$break_flag="user interrupt"   
                return()
              } else
              if (verbose_master) {print(paste0(Sys.time()," ...no breakfile detected.")); flush.console()}
            
              if (verbose_master) {print(paste0(Sys.time()," ...checking for execution timeout of slaves...")); flush.console()}
            
              if (!is.null(globvars$execution_timeout)) sleeptime=check_execution_timeout(verbose_master=verbose_master && difftime(Sys.time(), output_time, units="sec") > 10)
              if (verbose_master)
              {
                if (difftime(Sys.time(), output_time, units="sec") > 10)
                {
                  print(paste(Sys.time()," ...still waiting"));
                  flush.console()
                  output_time=Sys.time()
                }
              }
              
              if (verbose_master) {print(paste0(Sys.time()," ...going to sleep for ", sleeptime)); flush.console()}
              
              Sys.sleep(sleeptime) 
              
          }        
        if (!mpi.iprobe(mpi.any.source(),mpi.any.tag())) next #the previous loop was broken because of timeout, re-issue task
        if (verbose_master) {print(paste(Sys.time()," ...message detected")); flush.console()}
      slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        if (verbose_master) {print(paste(Sys.time()," ...message received")); flush.console()}
      slave_message_info <- mpi.get.sourcetag()
      slave_id <- slave_message_info[1]
      tag      <- slave_message_info[2]
        if (verbose_master) {print(paste(Sys.time()," ...message received from slave",slave_id, "(tag",tag,")")); flush.console()}          
#         browser()
     
      if (tag == 2) {      #retrieve result
          if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"sent results")); flush.console()}  
#        current_particle =which(globvars$node_id==slave_id & globvars$status==2)           #find which particle this result belongs to
         current_particle = globvars$slave_status[slave_id,"task_particle"]           #find which particle this result belongs to
  
        if (current_particle < 0)                                #redundant result of timeouted slave or its substitute
 				    {if (verbose_master) print(paste("slave",slave_id,"returned a now redundant result, ignored."))}
        else
        {
          globvars$fitness_X [current_particle] = slave_message
          globvars$status    [current_particle] =1      #mark as "finished"
          globvars$recent_error [current_particle] = 0  #mark as "no recent error"
          globvars$slave_status[slave_id,"timeouts_in_row"] = 0             #this result was in time, so reset the counter of consecutive timeouts 
          globvars$slave_status[slave_id,"task_particle"] = 0                     #remove association to completed particle
          redundant_calculations = globvars$slave_status[,"task_particle"] == current_particle  #find still treating the same particle
          if (any(redundant_calculations))
            globvars$slave_status[redundant_calculations, "task_particle"] = -current_particle #mark these slave as performing redundant calculations
          if (!is.null(globvars$execution_timeout)) 
            globvars$execution_times = rbind(globvars$execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),globvars$computation_start[current_particle],units="sec"))))   #monitor execution times

            if (init_search)
            {
              if (!is.null(globvars$execution_timeout)) globvars$execution_times = rbind(globvars$execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),globvars$computation_start[current_particle],units="sec"))))   #monitor execution times
              if (!is.null(logfile))  write.table(file = logfile, cbind(format(globvars$computation_start[current_particle],"%Y-%m-%d %H:%M:%S"), matrix(globvars$X[current_particle,],ncol=number_of_parameters)  , globvars$fitness_X[current_particle], 
              globvars$node_id[current_particle]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)          
            }
            else
            globvars$function_calls[current_particle] =globvars$function_calls[current_particle]+1        #increase iteration counter
          }
          if (globvars$slave_status[slave_id,"timeouts_in_row" ] < maxtries) #if this slave has not caused timeouts too often, keep on using it
            globvars$idle_slaves=c(globvars$idle_slaves,slave_id)
        } else if (tag == 3) {    # A slave has closed down.
            if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"has gone away")); flush.console()}
            globvars$closed_slaves = globvars$closed_slaves + 1
        } else if (tag == 4) {    # error occured during the execution of the objective function
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has produced an error:",slave_message))
            current_particle = globvars$slave_status[slave_id,"task_particle"]           #find which particle this result belongs to
            if (globvars$recent_error [current_particle] == 0)
              if (verbose_master) print(" ...retrying particle...") 
            else
            {   
              if (verbose_master) print(" This particle caused an error the second time, quitting...") 
              globvars$break_flag=paste("Abort, slave",slave_id,":",slave_message)
              globvars$closed_slaves=globvars$nslaves
            }  
            flush.console()
        } else if (tag == 5) {    #object request from slave 
    			slave_message_info <- mpi.get.sourcetag() 
    			 if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"requested object(s)",paste(slave_message, collapse=","))); flush.console()}
          obj_list=list()
          environment(get_object)=environment() 
          for (object_name in slave_message) #treat multiple requests, if required
             obj_list[[object_name]]= get_object(object_name=object_name) #get variable or parts thereof       

    			mpi.send.Robj(obj=obj_list, dest=slave_id, tag=5) #return object list to slave
    			 if (verbose_master) {print(paste(Sys.time()," ...object(s) sent to slave",slave_id)); flush.console()}
		    } else if (tag == 6) {    #object pushed from slave 
    			slave_message_info <- mpi.get.sourcetag() 
          if (globvars$slave_status[slave_id,"task_particle"] < 0) 
          {   #this is a now redundant call, ignore the push
               if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"is performing redundant call, push(es) ignored.")); flush.console()}
               next
          }  else
          environment(set_object)=environment() 
          for (i in 1:length(slave_message))
          {
      			object_name = names(slave_message[i]) #name of object
     			  if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"pushed object",object_name,"=",paste(slave_message[[i]], collapse=", "))); flush.console()}
            set_object(object_name=object_name, value=slave_message[[i]]) #set variable or parts thereof
          }
	      } else  {    #unknown tag
   			 if (verbose_master) {print(paste(Sys.time()," ...slave",slave_id,"sent unknown tag ",tag,", ignored")); flush.console()}
		    }  
  } 

  return()
}