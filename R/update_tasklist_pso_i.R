#internal function: update particle positions and velocities according to newly available results
update_tasklist_pso_i=function()                        #update particle positions for the next iteration based on all available results
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_* to allow accessing the same globals  (clumsy, but I don't know better)

{
   eval(parse(text=paste(c("do_plot_function=",     deparse(do_plot_function_i)))))  #this creates local version of the function do_plot_function 
   if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
      assign("break_flag","user interrupt",parent.frame())   
   
   completed_particles=status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no results available...don't update
   if (!all(completed_particles) & wait_complete_iteration) return()       #not all results available and complete iteration enabled...wait

    if (!is.null(logfile))        #append to logfile
      write.table(file = logfile, cbind(format(computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(X[completed_particles, ],ncol=ncol(X))  , fitness_X[completed_particles],
    node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

    do_plot_function()
   
    if (wait_for_keystroke && (!exists("ch") || ch!="c")) assign("ch",readline(),parent.frame()) 
 

   if (any(fitness_X %in% c(NA, NaN)))
      stop("Objective function mustn't yield NA nro NaN. Modify it to return very large numbers instead.")


   # Update the local bests and their fitness
   improved_particles=fitness_X < fitness_lbest #mark particles that improved their fitness
   if(is.na(completed_particles & improved_particles)[1]) #
   flush.console()

   if (any(improved_particles))
   {
     fitness_lbest[completed_particles & improved_particles] = fitness_X[completed_particles & improved_particles]            #store best fitness value
     X_lbest[completed_particles & improved_particles,1: number_of_parameters] = X[completed_particles & improved_particles,1: number_of_parameters]                  #store best parameter set
     assign("X_lbest",X_lbest,parent.frame())
     assign("fitness_lbest",fitness_lbest,parent.frame())
   }

#   if (!is.null(logfile))        #append to logfile
#    write.table(cbind(format(computation_start[completed_particles], "%Y-%m-%d %H:%M:%S"), X[completed_particles,], fitness_X[completed_particles],node_id[completed_particles]), file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

   # Update the global best and its fitness
   min_fitness_index = which.min(fitness_X[completed_particles])
   min_fitness =min(fitness_X[completed_particles])

   if (min_fitness < fitness_gbest)        #new global minimum found?
   {
       fitness_gbest = min_fitness          #update global minimum
       X_gbest[] = X[which(completed_particles)[min_fitness_index],]
       assign("X_gbest",X_gbest,parent.frame())
       assign("fitness_gbest",fitness_gbest,parent.frame())
       completed_particles = completed_particles | (status==0)     #also update the pending tasks with the new optimum
   }


   #check abortion criteria
   if ((fitness_itbest-fitness_gbest                   > abstol) &                     #check improvement in absolute terms
       abs((fitness_itbest-fitness_gbest)/max(0.00001,abs(fitness_gbest)) > reltol   ))  #check improvement in relative terms
   {    #improvement achieved
     it_last_improvent=max(function_calls)    #store iteration number which achieved the this improvement
     fitness_itbest=fitness_gbest
     assign("it_last_improvent",it_last_improvent,parent.frame())
     assign("fitness_itbest",fitness_itbest,parent.frame())
   } else
   {
      if (min(function_calls) - it_last_improvent>=max_wait_iterations)
        break_flag="converged"  #status=3
   }

   if (min(function_calls) >= max_number_of_iterations)
       break_flag="max iterations reached"

   if (!is.null(max_number_function_calls) && (sum(function_calls) >= max_number_function_calls))
       break_flag="max number of function calls reached"
   

    assign("evals_since_lastsave",evals_since_lastsave+sum(completed_particles),parent.frame())
    if ((evals_since_lastsave>=save_interval) | (!is.null(break_flag)))        #update project file if save interval is reached or if the abortion criterium is met
    {
      assign("evals_since_lastsave",0,parent.frame())
      if (!is.null(projectfile))
      {
		if (!is.null(colnames(X)))
		  par_names=colnames(X) else
		  par_names=paste(rep("par",number_of_parameters),seq(1,number_of_parameters),sep="_") #simple numbering of parameters

		  col.names=c(paste("best_",par_names,sep=""),"best_objective_function", paste("current_",par_names,sep=""),
          paste("current_velocity_",par_names,sep=""),"current_objective_function", "status", "begin_execution", "node_id","function_calls")
        write.table(cbind(X_lbest, fitness_lbest, X, V, fitness_X, status, format(computation_start, "%Y-%m-%d %H:%M:%S"), node_id, function_calls), file = projectfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
      }
      if(!is.null(plot_progress)) do.call(plot_optimization_progress, plot_progress)  #produce plots of optimization progress
    }

   if (!is.null(break_flag))
   {
      #assign("status",array(3,number_of_particles),parent.frame())   #set all particles to "done" -> abort optimization
      assign("break_flag",break_flag,parent.frame())   
      return()
   }

   
   # Update the particle velocity and position
   for (i in which(completed_particles))
   {
     R1 = runif(number_of_parameters)
     R2 = runif(number_of_parameters)
     V[i,] = w*V[i,] +
             C1*R1*(X_lbest[i,] - X[i,]) +
             C2*R2*(X_gbest - X[i,])
     V[i,] = V[i,] * min(1,abs(Vmax/V[i,]))        #limit to maximum velocity
     X[i,] = X[i,] + V[i,]
   }
   X[completed_particles, ] = t(pmax(t(X[completed_particles, ]),
        parameter_bounds[, 1]))
   X[completed_particles, ] = t(pmin(t(X[completed_particles, ]),
        parameter_bounds[, 2]))
   
   status   [completed_particles]=0      #mark as "to be computed"
   fitness_X[completed_particles]=Inf
   node_id  [completed_particles] =0

   #set globals
   assign("X",         X,parent.frame())
   assign("V",         V,parent.frame())
   assign("status",status,parent.frame())
   assign("node_id",node_id,parent.frame())
}
