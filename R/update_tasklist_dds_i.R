#internal function: update list of "particle positions" that need to be recorded
update_tasklist_dds <- function(loop_counter=1)                        
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{
   eval(parse(text=paste(c("do_plot_function=",     deparse(do_plot_function_i)))))  #this creates local version of the function do_plot_function 

   if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
      globvars$break_flag="user interrupt"  
   
   completed_particles=globvars$status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no new results available...don't update tasks

    if (!is.null(logfile) & loop_counter!=0)        #append to logfile, when enabled and when not in very first loop
      write.table(file = logfile, cbind(format(globvars$computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(globvars$X[completed_particles, ],ncol=ncol(globvars$X))  , globvars$fitness_X[completed_particles],
    globvars$node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

    do_plot_function()
   
    if (wait_for_keystroke && (!exists("globvars$ch") || globvars$ch!="c")) globvars$ch=readline()


   if (any(globvars$fitness_X %in% c(NA, NaN)))
      stop("Objective function mustn't yield NA nor NaN. Modify it to return very large numbers instead.")

   # Update the local bests and their fitness
   improved_particles=globvars$fitness_X < globvars$fitness_lbest #mark particles that improved their fitness

   globvars$futile_iter_count[!improved_particles]  = globvars$futile_iter_count[!improved_particles] + 1        #reset counter of futile iterations for improved particles
   globvars$futile_iter_count[ improved_particles]  = 0        #reset counter of futile iterations for improved particles     
   
   if (any(improved_particles))
   {
     globvars$fitness_lbest[completed_particles & improved_particles] = globvars$fitness_X[completed_particles & improved_particles]            #store best fitness value
     globvars$X_lbest[completed_particles & improved_particles,1: number_of_parameters] = globvars$X[completed_particles & improved_particles,1: number_of_parameters]                  #store best parameter set
   }
   
   # Update the global best and its fitness
   min_fitness_index = which.min(globvars$fitness_X[completed_particles])[1]
   min_fitness =min(globvars$fitness_X[completed_particles])

   if (min_fitness < globvars$fitness_gbest)        #new global minimum found?
   {
       globvars$fitness_gbest = min_fitness          #update global minimum
       globvars$X_gbest[] = globvars$X[which(completed_particles)[min_fitness_index],]
   }

   if (!exists('dds_ver')) dds_ver=2        #DDS-version (subtype for testing)

   if (number_of_particles > 1)
   {                                                                            #relocate "astray" particles
     if (dds_ver==1) toberelocated = (globvars$futile_iter_count==max(globvars$futile_iter_count)) & (globvars$fitness_lbest==max(globvars$fitness_lbest)) & (globvars$fitness_lbest!=Inf)    #1. find particles that are worst in both objective function AND improvement
     if (dds_ver==2) toberelocated = (globvars$fitness_lbest > globvars$fitness_gbest)     #2. relocate all but the best particle
     if (dds_ver==3) toberelocated = which.max(globvars$futile_iter_count * (globvars$fitness_lbest!=min(globvars$fitness_lbest)))[1]     #3. find particles is worst in improvement (but not the global best) and set to best improving particle

     if (any(toberelocated))
     {
       #cat(paste(sum(toberelocated),"particles globvars$relocated\n"))
       if (do_plot!=FALSE) globvars$relocated=cbind(globvars$X_lbest[toberelocated,],globvars$fitness_lbest[toberelocated])

       if ((dds_ver==1) || (dds_ver==2))  #version 1&2
       {
         globvars$X_lbest          [toberelocated,]= matrix(rep(globvars$X_gbest[],sum(toberelocated)),ncol=number_of_parameters, byrow = TRUE)         #relocate particles to current best
         globvars$fitness_lbest    [toberelocated] = globvars$fitness_gbest     #set to current best       #could also be set to particle with lowest globvars$futile_iter_count
         globvars$futile_iter_count[toberelocated] = 0             #zero iteration counter
       }
        if (dds_ver==3) #version 3
        {
         most_recent_improved = which.min(globvars$futile_iter_count)[1]       #get index to most recently improved particle
         globvars$X_lbest          [toberelocated,]= matrix(rep(globvars$X_lbest[most_recent_improved,],length(toberelocated)),ncol=number_of_parameters, byrow = TRUE)         #relocate particles to current best
         globvars$fitness_lbest    [toberelocated] = globvars$fitness_lbest[most_recent_improved]     #set to firness of particle with lowest globvars$futile_iter_count
         globvars$futile_iter_count[toberelocated] = globvars$futile_iter_count[most_recent_improved]             #copy iteration counter
        }

      }
   }

   #check abortion criteria
   if ((globvars$fitness_itbest-globvars$fitness_gbest                   > abstol) &                     #check improvement in absolute terms
       abs((globvars$fitness_itbest-globvars$fitness_gbest)/max(0.00001,abs(globvars$fitness_gbest)) > reltol   ))  #check improvement in relative terms
   {    #improvement achieved
     globvars$it_last_improvement=max(globvars$function_calls)    #store iteration number which achieved this improvement
     globvars$fitness_itbest=globvars$fitness_gbest
   } else
   {
      if (min(globvars$function_calls) - globvars$it_last_improvement>=max_wait_iterations)
        globvars$break_flag="converged"  #globvars$status=3
   }

#   if (min(globvars$function_calls) >= max_number_of_iterations)
#       globvars$break_flag="max globvars$function_calls reached"

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
    		par_names=paste(rep("par_",number_of_parameters),seq(1,number_of_parameters),sep="_") #simple numbering of parameters
        col.names=c(paste("best_",par_names,sep=""),"best_objective_function", paste("current_",par_names,sep=""),
          paste("current_velocity_",par_names,sep=""),"current_objective_function", "status", "begin_execution", "globvars$node_id","globvars$function_calls")
        write.table(file = projectfile, cbind(globvars$X_lbest, globvars$fitness_lbest, globvars$X, globvars$V, globvars$fitness_X, round(globvars$status), format(globvars$computation_start, "%Y-%m-%d %H:%M:%S"), globvars$node_id, globvars$function_calls + function_calls_init), quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
      }
      if(!is.null(plot_progress)) do.call(plot_optimization_progress, plot_progress)  #produce plots of optimization progress
    }

   if (!is.null(globvars$break_flag))
      return()

   # Update the particle position
   parameter_ranges=r*(parameter_bounds[,2]-parameter_bounds[,1]) #neighbourhood pertubation parameter * parameter range  
   for (i in which(completed_particles))
   {
#    p_inclusion=1-log(globvars$function_calls[i])/log(max_number_function_calls)    #probability of including a parameter in the search
    p_inclusion=1-log(globvars$function_calls[i])/log(max_number_function_calls/number_of_particles)    #probability of including a parameter in the search (parallel version)

    dimensions_to_search=NULL
    for (d in 1:number_of_parameters) 
    if (runif(1)<=p_inclusion)
      dimensions_to_search=c(dimensions_to_search,d)    #add parameter to "neighbourhood" to be searched
    if (is.null(dimensions_to_search)) dimensions_to_search=sample(number_of_parameters,1)     #include at least one parameter
    globvars$V[i,]=0

    globvars$V[i,dimensions_to_search] = rnorm(length(dimensions_to_search))*parameter_ranges[dimensions_to_search]
    globvars$X[i,] = globvars$X_lbest[i,] + globvars$V[i,]
    #reflect parameter if out-of-bounds
    params_below_bounds = globvars$X[i,] < parameter_bounds[,1]
    params_above_bounds = globvars$X[i,] > parameter_bounds[,2]

    {
      globvars$X[i,params_below_bounds]=parameter_bounds[params_below_bounds,1]+(parameter_bounds[params_below_bounds,1]-globvars$X[i,params_below_bounds])
      params_overreflected=params_below_bounds & (globvars$X[i,] > parameter_bounds[,2])    #identify parameters that are now above bounds due to the reflection
      globvars$X[i,params_overreflected]=parameter_bounds[params_overreflected,1]                                              #set to lower bounds
    }
    {
      globvars$X[i,params_above_bounds]=parameter_bounds[params_above_bounds,2]+(parameter_bounds[params_above_bounds,2]-globvars$X[i,params_above_bounds])
      params_overreflected=params_above_bounds & (globvars$X[i,] < parameter_bounds[,1])    #identify parameters that are now below bounds due to the reflection
      globvars$X[i,params_overreflected]=parameter_bounds[params_overreflected,2]                                              #set to upper bounds
    }
   }
   
   globvars$status   [completed_particles]=0      #mark as "to be computed"
   globvars$fitness_X[completed_particles]=Inf
   globvars$node_id  [completed_particles] =0

  if (!is.null(max_number_function_calls))
   {
		scheduled_calls = sum(globvars$status==0 | globvars$status==2)	#count how many tasks are scheduled or still running
        overcomitted_calls = (scheduled_calls +  sum(globvars$function_calls)) - max_number_function_calls #how many calls are scheduled that exceed max_number_function_calls
		if (overcomitted_calls > 0)
			globvars$status   [which(globvars$status==0)[1:overcomitted_calls]] = 0.1	#mark these as "ready to be calculated, but exceeding max_number_function_calls"
  }

}
