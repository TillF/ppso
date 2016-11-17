optim_pso <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40,max_number_of_iterations=5, max_number_function_calls=500, w=1,  C1=2, C2=2, abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   wait_complete_iteration=FALSE,parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3,lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="ppso.log",projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE,tryCall=FALSE, verbose=FALSE, ...)
# do particle swarm optimization
{
  
 # #algorithm parameters
#    number_of_particles=40
#    max_number_of_iterations=5
#    w=1                             #inertia constant
#    C1=2                            #cognitive components
#    C2=2                            #social component
#    #abort criteria
#    abstol=-Inf                      #minimum absolute improvement between iterations  (default: -Inf)
#    reltol=-Inf                       #minimum absolute improvement between iterations (default: -Inf)
#    max_wait_iterations=50    #number of iterations, within these an improvement of the above described quality has to be achieved (default:number_of_iterations)
#    #parallel options
#    wait_complete_iteration=FALSE    #wait for evaluation of all particles (complete iteration) before new iteration is started
#    objective_function=sample_function  #
#  
#  #problem parameters
#    number_of_parameters=2
#    parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters))     #matrix containing lower and uppaer boundary for each parameter
#    Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3    #maximum velocity
#  
#  #runtime & display parameters
#    do_plot=c(NULL,"base","rgl")             #enable 3D-plot of response surface and search progress (didactical purpose, only for two-parameter search and fast objective function. "base works only with  wait_complete_iteration=TRUE)
#    do_plot=NULL
#    wait_for_keystroke=FALSE                  #waiting for keystroke between iterations
#    logfile="log.txt"                         #logfile for optional logging of all model runs

#  

if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")
  
globvars$is_mpi = FALSE
  
environment(update_tasklist_pso)=environment() #force subroutines to have this function as parent (implicitly allowing read-only access to globals)
environment(init_particles)=environment() 
environment(init_visualisation)=environment() 
environment(plot_optimization_progress)=environment() 



if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   

globvars$evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected


#initialisation
globvars$nslaves = 1
globvars$X             =array(Inf,c(number_of_particles,number_of_parameters))  #globvars$X: position in parameter space                          
globvars$V             =array(0,c(number_of_particles,number_of_parameters))  #globvars$V: velocity in parameter space
globvars$fitness_X     =array(Inf,number_of_particles)            #fitness of each particle at current iteration
globvars$status        =array(0,number_of_particles)  #particle globvars$status: 0: to be computed; 1: finished; 2: in progress
globvars$computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if globvars$status=2)
globvars$node_id       =array(0,number_of_particles)                              #node number of worker / slave
globvars$X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
globvars$fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far
globvars$function_calls    =array(0,number_of_particles)  # iteration counter for each particle
                          
globvars$X_gbest     =array(Inf,number_of_parameters)            #global optimum

globvars$break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the local fitness to the worst possible
globvars$fitness_lbest[] = Inf
globvars$fitness_gbest = min(globvars$fitness_lbest);

init_particles(lhc_init)  #initialize velocities and particle positions

if (max_number_function_calls < 0)
{                                                         #indicator for "reset function counter" - ignore the number of function calls read from the project file
  globvars$function_calls[]=0
  max_number_function_calls=abs(max_number_function_calls)
}

if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  {
    if (!is.null(colnames(globvars$X)))
      par_names=paste(colnames(globvars$X),collapse="\t") else
      par_names=paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t") #simple numbering of parameters
    write.table(paste("time",par_names,"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }

 globvars$fitness_itbest= globvars$fitness_gbest     #best fitness in the last it_last iterations
 globvars$it_last_improvement=0               #counter for counting iterations since last improvement


while (is.null(globvars$break_flag))
{
  globvars$computation_start[] =Sys.time()      #store time of start of this computation loop
  if (wait_complete_iteration)      #evaluate all tasks before updating
  {    
      if (tryCall)                  #catch error message during evaluation (slower)
      {
        globvars$fitness_X=try(apply(globvars$X,1,objective_function, ...),silent=TRUE)
        if (!is.numeric(globvars$fitness_X))                      #an error occured during execution
        {
          globvars$break_flag=paste("aborted: ",as.character(globvars$fitness_X))    
          next
        }        
      }
      else
        globvars$fitness_X=apply(globvars$X,1,objective_function, ...)     #no error message during evaluation (faster)

      globvars$status    [] =1      #mark as "finished"
      globvars$node_id   [] =0
      globvars$function_calls[] =globvars$function_calls[]+1        #increase iteration counter
      update_tasklist_pso()   #update particle speeds and positions based on available results
  } else
  for (current_particle in 1:number_of_particles)      #do updates of tasks between single evaluations
  {
    if (tryCall)                  #catch error message during evaluation (slower)
    {
      globvars$fitness_X [current_particle]=try(objective_function(globvars$X[current_particle,], ...),silent=TRUE)
      if (!is.numeric(globvars$fitness_X [current_particle]))                      #an error occured during execution
      {
        globvars$break_flag=paste("aborted: ",as.character(globvars$fitness_X [current_particle]))    
        break
      }        
    }
    else
      globvars$fitness_X [current_particle] =objective_function(globvars$X[current_particle,], ...)     #no error message during evaluation (faster)

    globvars$status    [current_particle] =1      #mark as "finished"
    globvars$node_id   [current_particle] =0
    globvars$function_calls[current_particle] =globvars$function_calls[current_particle]+1        #increase iteration counter
    update_tasklist_pso()   #update particle speeds and positions based on available results
  }
    
}      
      
   
ret_val=list(value=globvars$fitness_gbest,par=globvars$X_gbest, function_calls=sum(globvars$function_calls),break_flag=globvars$break_flag) 

return(ret_val) 
}

