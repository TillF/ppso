#internal function: update list of "particle positions" that need to be recorded
update_tasklist_dds_i=function(loop_counter=1)                        
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_* to allow accessing the same globals  (clumsy, but I don't know better)

{
   if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
      assign("break_flag","user interrupt",parent.frame())   
   
   completed_particles=status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no new results available...don't update tasks

    if (!is.null(logfile) & loop_counter!=0)        #append to logfile, when enabled and when not in very first loop
      write.table(file = logfile, cbind(format(computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(X[completed_particles, ],ncol=ncol(X))  , fitness_X[completed_particles],
    node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)


  if (("base" %in% do_plot))       #do base plotting
    {
      if (exists("plot_window") && (plot_window %in% dev.list()))     #activate progress_plot window, if already open
      {
        dev.set(plot_window)
      }  else
      {
        windows()
        assign("plot_window", dev.cur(), pos=parent.frame())
      }

      res = persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",     ltheta = 120, shade = 0.75, ticktype = "detailed",      xlab = "X", ylab = "Y", zlab = "obj fun")
      points(trans3d(X[,1], X[,2], fitness_X, pmat = res), col = 2, pch =16)
    }

    if ("rgl" %in% do_plot)        #do rgl plotting
    {
      rgl.pop(id=hdl[completed_particles & hdl!=0])     #remove outdated dots
      for (i in which(completed_particles))
         hdl[i]=points3d(X[i,1], X[i,2], fitness_X[i], col="red")
      assign("hdl",hdl,parent.frame())
    }

    if (wait_for_keystroke) readline()


   # Update the local bests and their fitness

   improved_particles=fitness_X < fitness_lbest #mark particles that improved their fitness

   if (any(improved_particles))
   {
     fitness_lbest[completed_particles & improved_particles] = fitness_X[completed_particles & improved_particles]            #store best fitness value
     X_lbest[completed_particles & improved_particles,1: number_of_parameters] = X[completed_particles & improved_particles,1: number_of_parameters]                  #store best parameter set
     assign("X_lbest",X_lbest,parent.frame())
     assign("fitness_lbest",fitness_lbest,parent.frame())
   }

   # Update the global best and its fitness
   min_fitness_index = which.min(fitness_X[completed_particles])
   min_fitness =min(fitness_X[completed_particles])

   if (min_fitness < fitness_gbest)        #new global minimum found?
   {
       fitness_gbest = min_fitness          #update global minimum
       X_gbest[] = X[which(completed_particles)[min_fitness_index],]
       assign("X_gbest",X_gbest,parent.frame())
       assign("fitness_gbest",fitness_gbest,parent.frame())
   }

   #check abortion criteria
   if ((fitness_itbest-fitness_gbest                   > abstol) &                     #check improvement in absolute terms
       abs((fitness_itbest-fitness_gbest)/max(0.00001,abs(fitness_gbest)) > reltol   ))  #check improvement in relative terms
   {    #improvement achieved
     it_last_improvent=max(iterations)    #store iteration number which achieved this improvement
     fitness_itbest=fitness_gbest
     assign("it_last_improvent",it_last_improvent,parent.frame())
     assign("fitness_itbest",fitness_itbest,parent.frame())
   } else
   {
      if (min(iterations) - it_last_improvent>=max_wait_iterations)
        break_flag="converged"  #status=3
   }

#   if (min(iterations) >= max_number_of_iterations)
#       break_flag="max iterations reached"

   if (!is.null(max_number_function_calls) && (sum(iterations) >= max_number_function_calls))
       break_flag="max number of function calls reached"

  
    assign("evals_since_lastsave",evals_since_lastsave+sum(completed_particles),parent.frame())
    if ((evals_since_lastsave>=save_interval) | (!is.null(break_flag)))        #update project file if save interval is reached or if the abortion criterium is met
    {
      assign("evals_since_lastsave",0,parent.frame())
      if (!is.null(projectfile))
      {
        col.names=c(paste("best_par_",1:ncol(X),sep=""),"best_objective_function", paste("current_par_",1:ncol(X),sep=""),
          paste("current_velocity_par_",1:ncol(X),sep=""),"current_objective_function", "status", "begin_execution", "node_id","function_calls")
        write.table(cbind(X_lbest, fitness_lbest, X, V, fitness_X, status, format(computation_start, "%Y-%m-%d %H:%M:%S"), node_id, iterations), file = projectfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
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
#    p_inclusion=1-log(iterations[i])/log(max_number_function_calls)    #probability of including a parameter in the search
    p_inclusion=1-log(iterations[i])/log(max_number_function_calls/number_of_particles)    #probability of including a parameter in the search (parallel version)

    dimensions_to_search=NULL
    for (d in 1:number_of_parameters) 
    if (runif(1)<=p_inclusion)
      dimensions_to_search=c(dimensions_to_search,d)    #add parameter to "neighbourhood" to be searched
    if (is.null(dimensions_to_search)) dimensions_to_search=sample(number_of_parameters,1)     #include at least one parameter
    V[i,]=0
#    browser()
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

   #set globals
   assign("X",         X,parent.frame())
   assign("V",         V,parent.frame())
   assign("status",status,parent.frame())
   assign("node_id",node_id,parent.frame())
}
