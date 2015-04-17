#internal function: update particle positions and velocities according to newly available results
update_tasklist_pso=function()                        #update particle positions for the next iteration based on all available results
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{
   environment(do_plot_function)=environment()  #this creates local version of the function do_plot_function 
   if (!is.null(max_number_function_calls) && (sum(globvars$function_calls) >= max_number_function_calls))
       globvars$break_flag="max number of function calls reached"

   if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
      globvars$break_flag="user interrupt"  
   
   completed_particles=globvars$status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no results available...don't update
   if (!all(completed_particles) & wait_complete_iteration) return()       #not all results available and complete iteration enabled...wait

    if (!is.null(logfile))        #append to logfile
      write.table(file = logfile, cbind(format(globvars$computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(globvars$X[completed_particles, ],ncol=ncol(globvars$X))  , globvars$fitness_X[completed_particles],
    globvars$node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

    do_plot_function()
   
    if (wait_for_keystroke && (!exists("ch", where=globvars) || globvars$ch!="c")) 
    {
      if (is.numeric(globvars$ch) && (globvars$ch > 0)) globvars$ch=globvars$ch-1 else
      {
        print("ENTER to proceed, 'b'+ENTER for debug mode, <n> to skip <n> times, 'c'+ENTER to continue till end")
        globvars$ch=readline() 
      }  
      if (globvars$ch=="b")	browser() else
      if (globvars$ch!="c") globvars$ch = strtoi(globvars$ch) #counter
      if (is.na(globvars$ch)) globvars$ch=""
    }  
   if (any(globvars$fitness_X %in% c(NA, NaN)))
      stop("Objective function mustn't yield NA or NaN. Modify it to return very large numbers instead.")


   # Update the local bests and their fitness
   improved_particles=globvars$fitness_X < globvars$fitness_lbest #mark particles that improved their fitness
   if(is.na(completed_particles & improved_particles)[1]) #
   flush.console()

   if (any(improved_particles))
   {
     globvars$fitness_lbest[completed_particles & improved_particles] = globvars$fitness_X[completed_particles & improved_particles]            #store best fitness value
     globvars$X_lbest[completed_particles & improved_particles,1: number_of_parameters] = globvars$X[completed_particles & improved_particles,1: number_of_parameters]                  #store best parameter set
   }

#   if (!is.null(logfile))        #append to logfile
#    write.table(cbind(format(globvars$computation_start[completed_particles], "%Y-%m-%d %H:%M:%S"), globvars$X[completed_particles,], globvars$fitness_X[completed_particles],globvars$node_id[completed_particles]), file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

   # Update the global best and its fitness
   min_fitness_index = which.min(globvars$fitness_X[completed_particles])
   min_fitness =min(globvars$fitness_X[completed_particles])

   if (min_fitness < globvars$fitness_gbest)        #new global minimum found?
   {
       globvars$fitness_gbest = min_fitness          #update global minimum value
       globvars$X_gbest[] = globvars$X[which(completed_particles)[min_fitness_index],] #update global minimum
       completed_particles = completed_particles | ((globvars$status==0) & (globvars$node_id != -1))     #also update the pending tasks (which are not fixed (initial estimates)) with the new optimum
   }


   #check abortion criteria
   if ((globvars$fitness_itbest-globvars$fitness_gbest                   > abstol) &                     #check improvement in absolute terms
       abs((globvars$fitness_itbest-globvars$fitness_gbest)/max(0.00001,abs(globvars$fitness_gbest)) > reltol   ))  #check improvement in relative terms
   {    #improvement achieved
     globvars$it_last_improvement=max(globvars$function_calls)    #store iteration number which achieved the this improvement
     globvars$fitness_itbest=globvars$fitness_gbest
   } else
   {
      if (min(globvars$function_calls) - globvars$it_last_improvement>=max_wait_iterations)
        globvars$break_flag="converged"  #globvars$status=3
   }

   if (min(globvars$function_calls) >= max_number_of_iterations)
       globvars$break_flag="max iterations reached"

   if (!is.null(max_number_function_calls) && (sum(globvars$function_calls) >= max_number_function_calls))
       globvars$break_flag="max number of function calls reached"
   

    globvars$evals_since_lastsave = globvars$evals_since_lastsave+sum(completed_particles)
    if ((globvars$evals_since_lastsave>=save_interval) | (!is.null(globvars$break_flag)))        #update project file if save interval is reached or if the abortion criterium is met
    {
      globvars$evals_since_lastsave=0
      if (!is.null(projectfile))
      {
		if (!is.null(colnames(globvars$X)))
		  par_names=colnames(globvars$X) else
		  par_names=paste(rep("par",number_of_parameters),seq(1,number_of_parameters),sep="_") #simple numbering of parameters

		  col.names=c(paste("best_",par_names,sep=""),"best_objective_function", paste("current_",par_names,sep=""),
          paste("current_velocity_",par_names,sep=""),"current_objective_function", "status", "begin_execution", "node_id","function_calls")
        write.table(cbind(globvars$X_lbest, globvars$fitness_lbest, globvars$X, globvars$V, globvars$fitness_X,globvars$status, format(globvars$computation_start, "%Y-%m-%d %H:%M:%S"), globvars$node_id, globvars$function_calls),
          file = projectfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
      }
      if(!is.null(plot_progress)) do.call(plot_optimization_progress, plot_progress)  #produce plots of optimization progress
    }

   if (!is.null(globvars$break_flag))
      return()

   
   # Update the particle velocity and position
  globvars$status   [completed_particles]=0      #mark as "to be computed"
  globvars$fitness_X[completed_particles]=Inf
  globvars$node_id  [completed_particles] =0

    for (i in which(completed_particles))
   {
     if (ncol(globvars$pending_initial_estimates) > 0) #there are still initial estimates that need to be evaluated
     {
       globvars$X[i,] = globvars$pending_initial_estimates[,1]
       globvars$pending_initial_estimates = globvars$pending_initial_estimates[,-1, drop=FALSE] #remove from list
       globvars$node_id   [i]= -1      #mark as "dont change anymore until finished"
     } else
     {
       R1 = runif(number_of_parameters)
       R2 = runif(number_of_parameters)
       globvars$V[i,] = w*globvars$V[i,] +
               C1*R1*(globvars$X_lbest[i,] - globvars$X[i,]) +
               C2*R2*(globvars$X_gbest - globvars$X[i,])
       globvars$V[i,] = globvars$V[i,] * min(1,abs(Vmax/globvars$V[i,]))        #limit to maximum velocity
       globvars$X[i,] = globvars$X[i,] + globvars$V[i,]
     }
   }
   globvars$X[completed_particles, ] = t(pmax(t(globvars$X[completed_particles, ]),
        parameter_bounds[, 1]))
   globvars$X[completed_particles, ] = t(pmin(t(globvars$X[completed_particles, ]),
        parameter_bounds[, 2]))
   
   

}
