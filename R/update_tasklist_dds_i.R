#internal function: update list of "particle positions" that need to be recorded
update_tasklist_dds_i <- function(loop_counter=1)                        
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_* to allow accessing the same globals  (clumsy, but I don't know better)

{
   eval(parse(text=paste(c("do_plot_function=",     deparse(do_plot_function_i)))))  #this creates local version of the function do_plot_function 

   if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
      assign("break_flag","user interrupt",parent.frame())   
   
   completed_particles=status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no new results available...don't update tasks

    if (!is.null(logfile) & loop_counter!=0)        #append to logfile, when enabled and when not in very first loop
      write.table(file = logfile, cbind(format(computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(X[completed_particles, ],ncol=ncol(X))  , fitness_X[completed_particles],
    node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

    do_plot_function()
   
    if (wait_for_keystroke && (!exists("ch") || ch!="c")) assign("ch",readline(),parent.frame()) 


   # Update the local bests and their fitness

   improved_particles=fitness_X < fitness_lbest #mark particles that improved their fitness

   futile_iter_count[!improved_particles]  = futile_iter_count[!improved_particles] + 1        #reset counter of futile iterations for improved particles
   futile_iter_count[ improved_particles]  = 0        #reset counter of futile iterations for improved particles     
   assign("futile_iter_count",futile_iter_count,parent.frame())
   
   if (any(improved_particles))
   {
     fitness_lbest[completed_particles & improved_particles] = fitness_X[completed_particles & improved_particles]            #store best fitness value
     X_lbest[completed_particles & improved_particles,1: number_of_parameters] = X[completed_particles & improved_particles,1: number_of_parameters]                  #store best parameter set
     assign("X_lbest",X_lbest,parent.frame())
     assign("fitness_lbest",fitness_lbest,parent.frame())
   }
   
   # Update the global best and its fitness
   min_fitness_index = which.min(fitness_X[completed_particles])[1]
   min_fitness =min(fitness_X[completed_particles])

   if (min_fitness < fitness_gbest)        #new global minimum found?
   {
       fitness_gbest = min_fitness          #update global minimum
       X_gbest[] = X[which(completed_particles)[min_fitness_index],]
       assign("X_gbest",X_gbest,parent.frame())
       assign("fitness_gbest",fitness_gbest,parent.frame())
   }

   if (!exists('dds_ver')) dds_ver=2        #DDS-Version (subtype for testing)

   if (number_of_particles > 1)
   {                                                                            #relocate "astray" particles
     if (dds_ver==1) toberelocated = (futile_iter_count==max(futile_iter_count)) & (fitness_lbest==max(fitness_lbest)) & (fitness_lbest!=Inf)    #1. find particles that are worst in both objective function AND improvement
     if (dds_ver==2) toberelocated = (fitness_lbest > fitness_gbest)     #2. relocate all but the best particle
     if (dds_ver==3) toberelocated = which.max(futile_iter_count * (fitness_lbest!=min(fitness_lbest)))[1]     #3. find particles is worst in improvement (but not the global best) and set to best improving particle

     if (any(toberelocated))
     {
       #cat(paste(sum(toberelocated),"particles relocated\n"))
       if (do_plot!=FALSE) assign("relocated",cbind(X_lbest[toberelocated,],fitness_lbest[toberelocated]),parent.frame())

       if ((dds_ver==1) || (dds_ver==2))  #version 1&2
       {
         X_lbest          [toberelocated,]= matrix(rep(X_gbest[],sum(toberelocated)),ncol=number_of_parameters, byrow = TRUE)         #relocate particles to current best
         fitness_lbest    [toberelocated] = fitness_gbest     #set to current best       #could also be set to particle with lowest futile_iter_count
         futile_iter_count[toberelocated] = 0             #zero iteration counter
       }
        if (dds_ver==3) #version 3
        {
         most_recent_improved = which.min(futile_iter_count)[1]       #get index to most recently improved particle
         X_lbest          [toberelocated,]= matrix(rep(X_lbest[most_recent_improved,],length(toberelocated)),ncol=number_of_parameters, byrow = TRUE)         #relocate particles to current best
         fitness_lbest    [toberelocated] = fitness_lbest[most_recent_improved]     #set to firness of particle with lowest futile_iter_count
         futile_iter_count[toberelocated] = futile_iter_count[most_recent_improved]             #copy iteration counter
        }

       assign("X_lbest",X_lbest,parent.frame())       
       assign("fitness_lbest",fitness_lbest,parent.frame())
       assign("futile_iter_count",futile_iter_count,parent.frame())
      }
   }

   #check abortion criteria
   if ((fitness_itbest-fitness_gbest                   > abstol) &                     #check improvement in absolute terms
       abs((fitness_itbest-fitness_gbest)/max(0.00001,abs(fitness_gbest)) > reltol   ))  #check improvement in relative terms
   {    #improvement achieved
     it_last_improvement=max(function_calls)    #store iteration number which achieved this improvement
     fitness_itbest=fitness_gbest
     assign("it_last_improvement",it_last_improvement,parent.frame())
     assign("fitness_itbest",fitness_itbest,parent.frame())
   } else
   {
      if (min(function_calls) - it_last_improvement>=max_wait_iterations)
        break_flag="converged"  #status=3
   }

#   if (min(function_calls) >= max_number_of_iterations)
#       break_flag="max function_calls reached"

   if (!is.null(max_number_function_calls) && (sum(function_calls) >= max_number_function_calls))
       break_flag="max number of function calls reached"

  
    assign("evals_since_lastsave",evals_since_lastsave+sum(completed_particles),parent.frame())
    if ((evals_since_lastsave>=save_interval) | (!is.null(break_flag)))        #update project file if save interval is reached or if the abortion criterium is met
    {
      assign("evals_since_lastsave",0,parent.frame())
      if (!is.null(projectfile))
      {
        col.names=c(paste("best_par_",1:ncol(X),sep=""),"best_objective_function", paste("current_par_",1:ncol(X),sep=""),
          paste("current_velocity_par_",1:ncol(X),sep=""),"current_objective_function", "status", "begin_execution", "node_id","function_calls")
        write.table(file = projectfile, cbind(X_lbest, fitness_lbest, X, V, fitness_X, round(status), format(computation_start, "%Y-%m-%d %H:%M:%S"), node_id, function_calls + function_calls_init), quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
      }
      if(!is.null(plot_progress)) do.call(plot_optimization_progress, plot_progress)  #produce plots of optimization progress
    }

   if (!is.null(break_flag))
   {
      assign("break_flag",break_flag,parent.frame())   
      return()
   }

   # Update the particle position
   parameter_ranges=r*(parameter_bounds[,2]-parameter_bounds[,1]) #neighbourhood pertubation parameter * parameter range  
   for (i in which(completed_particles))
   {
#    p_inclusion=1-log(function_calls[i])/log(max_number_function_calls)    #probability of including a parameter in the search
    p_inclusion=1-log(function_calls[i])/log(max_number_function_calls/number_of_particles)    #probability of including a parameter in the search (parallel version)

    dimensions_to_search=NULL
    for (d in 1:number_of_parameters) 
    if (runif(1)<=p_inclusion)
      dimensions_to_search=c(dimensions_to_search,d)    #add parameter to "neighbourhood" to be searched
    if (is.null(dimensions_to_search)) dimensions_to_search=sample(number_of_parameters,1)     #include at least one parameter
    V[i,]=0

    V[i,dimensions_to_search] = rnorm(length(dimensions_to_search))*parameter_ranges[dimensions_to_search]
    X[i,] = X_lbest[i,] + V[i,]
    #reflect parameter if out-of-bounds
    params_below_bounds = X[i,] < parameter_bounds[,1]
    params_above_bounds = X[i,] > parameter_bounds[,2]

    {
      X[i,params_below_bounds]=parameter_bounds[params_below_bounds,1]+(parameter_bounds[params_below_bounds,1]-X[i,params_below_bounds])
      params_overreflected=params_below_bounds & (X[i,] > parameter_bounds[,2])    #identify parameters that are now above bounds due to the reflection
      X[i,params_overreflected]=parameter_bounds[params_overreflected,1]                                              #set to lower bounds
    }
    {
      X[i,params_above_bounds]=parameter_bounds[params_above_bounds,2]+(parameter_bounds[params_above_bounds,2]-X[i,params_above_bounds])
      params_overreflected=params_above_bounds & (X[i,] < parameter_bounds[,1])    #identify parameters that are now below bounds due to the reflection
      X[i,params_overreflected]=parameter_bounds[params_overreflected,2]                                              #set to upper bounds
    }
   }
   
   status   [completed_particles]=0      #mark as "to be computed"
   fitness_X[completed_particles]=Inf
   node_id  [completed_particles] =0

  if (!is.null(max_number_function_calls)
   {
		scheduled_calls = sum(status==0 | status==2)	#count how many tasks are scheduled or still running
        overcomitted_calls = max_number_function_calls - (scheduled_calls +  sum(function_calls)) #how many calls are scheduled that exceed max_number_function_calls
		if (overcomitted_calls > 0)
			status   [which(status==0)[1:overcomitted_calls]] = 0.1	#mark these as "ready to be calculated, but exceeding max_number_function_calls"
  }

   

   #set globals
   assign("X",         X,parent.frame())
   assign("V",         V,parent.frame())
   assign("status",status,parent.frame())
   assign("node_id",node_id,parent.frame())
}
