#internal function: update particle positions and velocities according to newly available results
update_tasklist_pso_i=function()                        #update particle positions for the next iteration based on all available results
#note: this function reads and writes to non-local variables (i.e. varaibles declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_p*  (clumsy, but I don't know better)

{
   if ((!is.null(break_file)) & (file.exists(break_file)))      #check if interrupt by user is requested
      assign("break_flag","user interrupt",parent.frame())   
   
   completed_particles=status==1                   #mark completed particles
   if (all(completed_particles==FALSE)) return()                   #no results available...don't update
   if (!all(completed_particles) & wait_complete_iteration) return()       #not all results available and complete iteration enabled...wait

    if (!is.null(logfile))        #append to logfile
      write.table(file = logfile, cbind(format(computation_start[completed_particles],"%Y-%m-%d %H:%M:%S"), matrix(X[completed_particles, ],ncol=ncol(X))  , fitness_X[completed_particles],
    node_id[completed_particles]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)


  if (("base" %in% do_plot) & wait_complete_iteration)       #do base plotting
    {
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
     it_last_improvent=max(iterations)    #store iteration number which achieved the this improvement
     fitness_itbest=fitness_gbest
     assign("it_last_improvent",it_last_improvent,parent.frame())
     assign("fitness_itbest",fitness_itbest,parent.frame())
   } else
   {
      if (min(iterations) - it_last_improvent>=max_wait_iterations)
        break_flag="converged"  #status=3
   }

   if (min(iterations) >= max_number_of_iterations)
       break_flag="max iterations reached"

   

    assign("evals_since_lastsave",evals_since_lastsave+sum(completed_particles),parent.frame())
    if ((evals_since_lastsave>=save_interval) | (!is.null(break_flag)))        #update project file if save interval is reached or if the abortion criterium is met
    {
      assign("evals_since_lastsave",0,parent.frame())
      if (!is.null(projectfile))
      {
        col.names=c(paste("best_par_",1:ncol(X),sep=""),"best_objective_function", paste("current_par_",1:ncol(X),sep=""),
        paste("current_velocity_par_",1:ncol(X),sep=""),"current_objective_function", "status", "begin_execution", "node_id")
        write.table(cbind(X_lbest, fitness_lbest, X, V, fitness_X, status, format(computation_start, "%Y-%m-%d %H:%M:%S"), node_id), file = projectfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)
      }
    }

   if (!is.null(break_flag))
   {
      #assign("status",array(3,number_of_particles),parent.frame())   #set all particles to "done" -> abort optimization
      assign("break_flag",break_flag,parent.frame())   
      return()
   }

   if(max(which(completed_particles))>2) browser()
   
   # Update the particle velocity and position
   for (i in which(completed_particles))
   {
     R1 = runif(number_of_parameters)
     R2 = runif(number_of_parameters)
     V[i,] = w*V[i,] +
             C1*R1*(X_lbest[i,] - X[i,]) +
             C2*R2*(X_gbest - X[i,])
     if (all(V==0)) V[i,]= runif(min=-0.01, max=0.01)*Vmax else       #ensure that a particle doesn't stand still
     V[i,] = V[i,] * min(1,abs(Vmax/V[i,]))        #limit to maximum velocity
     X[i,] = X[i,] + V[i,]
   }
   X[completed_particles,]=pmax(X[completed_particles,],parameter_bounds[,1])    #restrain by specified boundaries
   X[completed_particles,]=pmin(X[completed_particles,],parameter_bounds[,2])
   
   status   [completed_particles]=0      #mark as "to be computed"
   fitness_X[completed_particles]=Inf
   node_id  [completed_particles] =0

   #set globals
   assign("X",         X,parent.frame())
   assign("V",         V,parent.frame())
   assign("status",status,parent.frame())
   assign("node_id",node_id,parent.frame())
}
