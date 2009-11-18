optim_pso <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40,max_number_of_iterations=5, max_number_function_calls=NULL, w=1,  C1=2, C2=2, abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   wait_complete_iteration=FALSE,parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3,lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="ppso.log",projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE,tryCall=FALSE)
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

  
eval(parse(text=paste(c("update_tasklist_pso=",deparse(update_tasklist_pso_i)))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles=",     deparse(init_particles_i)))))  #this creates local version of the function init_particles (see explanation there)
eval(parse(text=paste(c("init_visualisation=",     deparse(init_visualisation_i)))))  #this creates local version of the function init_visualisation 
if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   


evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected


#initialisation
X             =array(0,c(number_of_particles,number_of_parameters))  #X: position in parameter space                          
V             =array(0,c(number_of_particles,number_of_parameters))  #V: velocity in parameter space
fitness_X     =array(Inf,number_of_particles)            #fitness of each particle at current iteration
status        =array(0,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
node_id       =array(0,number_of_particles)                              #node number of worker / slave
X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far
iterations    =array(0,number_of_particles)  # iteration counter for each particle
                          
X_gbest     =array(Inf,number_of_parameters)            #global optimum

break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the local fitness to the worst possible
fitness_lbest[] = Inf
fitness_gbest = min(fitness_lbest);

init_particles(lhc_init)  #initialize velocities and particle positions



if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  write.table(paste("time",paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t"),"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


 fitness_itbest= fitness_gbest     #best fitness in the last it_last iterations
 it_last_improvent=0               #counter for counting iterations since last improvement


while (is.null(break_flag))
{
  if (wait_complete_iteration)      #evaluate all tasks before updateing
  {    
      if (tryCall)                  #catch error message during evaluation (slower)
      {
        fitness_X=try(apply(X,1,objective_function),silent=TRUE)
        if (!is.numeric(fitness_X))                      #an error occured during execution
        {
          break_flag=paste("aborted: ",as.character(fitness_X))    
          next
        }        
      }
      else
        fitness_X=apply(X,1,objective_function)     #no error message during evaluation (faster)

      status    [] =1      #mark as "finished"
      iterations[] =iterations[]+1        #increase iteration counter
      update_tasklist_pso()   #update particle speeds and positions based on available results
  } else
  for (current_particle in 1:number_of_particles)      #do updates of tasks between single evaluations
  {
    if (tryCall)                  #catch error message during evaluation (slower)
    {
      fitness_X [current_particle]=try(objective_function(X[current_particle,]),silent=TRUE)
      if (!is.numeric(fitness_X [current_particle]))                      #an error occured during execution
      {
        break_flag=paste("aborted: ",as.character(fitness_X [current_particle]))    
        break
      }        
    }
    else
      fitness_X [current_particle] =objective_function(X[current_particle,])     #no error message during evaluation (faster)
  
    status    [current_particle] =1      #mark as "finished"
    iterations[current_particle] =iterations[current_particle]+1        #increase iteration counter
    update_tasklist_pso()   #update particle speeds and positions based on available results
  }
    
}      
      
   
ret_val=list(value=fitness_gbest,par=X_gbest,iterations=min(iterations),break_flag=break_flag) 

return(ret_val) 
}

