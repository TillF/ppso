optim_dds <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles= 1, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)),lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, tryCall=FALSE)
# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{
  
 # #algorithm parameters
#    number_of_particles=40         #number of DDS-thread that are tracked
#    r=0.2                             #neighbourhood size perturbation parameter 
#    #abort criteria
#    abstol=-Inf                      #minimum absolute improvement between iterations  (default: -Inf)
#    reltol=-Inf                       #minimum absolute improvement between iterations (default: -Inf)
#    max_wait_iterations=50    #number of iterations, within these an improvement of the above described quality has to be achieved (default:number_of_iterations)
#    #parallel options
#    objective_function=sample_function  #
#  
#  #problem parameters
#    number_of_parameters=2
#    parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters))     #matrix containing lower and uppaer boundary for each parameter
#  
#  #runtime & display parameters
#    do_plot=c(NULL,"base","rgl")             #enable 3D-plot of response surface and search progress (didactical purpose, only for two-parameter search and fast objective function. "base works only with  wait_complete_iteration=TRUE)
#    do_plot=NULL
#    wait_for_keystroke=FALSE                  #waiting for keystroke between iterations
#    logfile="log.txt"                         #logfile for optional logging of all model runs

#  

if (!is.null(max_number_function_calls) && max_number_function_calls < number_of_particles)
  stop("max_number_function_calls must be at least number_of_particles.")
  
eval(parse(text=paste(c("update_tasklist_dds=",deparse(update_tasklist_dds_i)))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles=",     deparse(init_particles_i)))))  #this creates local version of the function init_particles (see explanation there)
eval(parse(text=paste(c("init_visualisation=",     deparse(init_visualisation_i)))))  #this creates local version of the function init_visualisation 

if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   


evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected

#initialisation
init_calls=ceiling(max(0.005*max_number_function_calls,5)) #number of function calls used to initialise each particle, if no value can be loaded from a project file
number_of_particles_org=number_of_particles                 #save original number of particles
number_of_particles=number_of_particles*init_calls          #increase number of particles for pre-search


nslaves = 1
X             =array(Inf,c(number_of_particles,number_of_parameters))  #X: position in parameter space                          
V             =array(0,c(           number_of_particles,number_of_parameters))  #V: velocity in parameter space
fitness_X     =array(Inf,number_of_particles)            #optimum of each particle at current iteration
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

#presearch / initialisation: 
#  the particles are preferrably initialized with the data from the projectfile. If that does not exist or does not contain enough records,
#  for each uninitialized particle (uninitialized_particles) a number of prior calls (init_calls) are performed, of which the best is used
  init_particles(lhc_init)  #initialize particle positions
  if (max_number_function_calls < 0)
  {                                                         #indicator for "reset function counter" - ignore the number of function calls read from the project file
    iterations[]=0
    max_number_function_calls=abs(max_number_function_calls)
  }
  
  if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
    write.table(paste("time",paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t"),"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  #restore array dimensions according to original number of particles
  number_of_particles=number_of_particles_org         #back to original number of particles
  X_lbest       =matrix(X_lbest      [1:number_of_particles,],ncol=number_of_parameters)        # current optimum of each particle so far
  fitness_lbest =       fitness_lbest[1:number_of_particles]                                    #best solution for each particle so far
  
  offseti=(0:(number_of_particles-1))*init_calls      #for extracting every init_calls'th element from the pre-run results (which contain up to number_of_particles*init_calls runs)
  uninitialized_particles= which(fitness_lbest==Inf)                      #"real" particles that still need to be initialized with a function value
  pre_run_computations = rep(1:number_of_particles,init_calls) %in% uninitialized_particles  #"trial" particles that need their function value to be computed (this results in flagging "init_calls" runs for every particle not yet initialized)
  if (sum(pre_run_computations) >= max_number_function_calls) stop(paste("Parameter max_number_function_calls =",max_number_function_calls,"does not suffice for initialisation. Increase it or decrease number_of_particles"))
 
  if (any(pre_run_computations))
  {
    computation_start[pre_run_computations]=Sys.time()
    fitness_X [pre_run_computations]=apply(X[pre_run_computations,],1,objective_function)    #execute pending pre-runs
    status    [pre_run_computations] =1      #mark as pre-runs "finished"
  
    if (!is.null(logfile))  write.table(file = logfile, cbind(format(computation_start[pre_run_computations],"%Y-%m-%d %H:%M:%S"), matrix(X[pre_run_computations,],ncol=ncol(X)), fitness_X[pre_run_computations], #write pre-runs to logfile, too
    node_id[pre_run_computations]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)
  } 

  max_number_function_calls=max_number_function_calls-sum(pre_run_computations)  #reduce number of available calls due to pre-search


  for (i in uninitialized_particles)  #initialize each uninitialized particle with the best of its "init_calls" pre-runs
  {
    min_fitness_index = which.min(fitness_X[i+offseti])            
    fitness_lbest[i] =fitness_X [i+offseti[min_fitness_index]]
    X_lbest      [i,]=X         [i+offseti[min_fitness_index],]
  }
  
  #restore array dimensions according to original number of particles
  X             =X_lbest  #X: position in parameter space                          
  V             =matrix(V[1:number_of_particles,],ncol=number_of_parameters)   #V: velocity in parameter space
  fitness_X     =fitness_X[1:number_of_particles]            #optimum of each particle at current iteration
#  status        =array(1,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
  computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
  node_id       =array(0,number_of_particles)                              #node number of worker / slave
  iterations    =iterations[1:number_of_particles]  # iteration counter for each particle
  status=status_org[1:number_of_particles]  #  restore original contents 
  futile_iter_count = futile_iter_count[1:number_of_particles]
  fitness_gbest = min(fitness_lbest);

# actual search

fitness_itbest= fitness_gbest     #best fitness in the last it_last iterations
it_last_improvent=0               #counter for counting iterations since last improvement

update_tasklist_dds(loop_counter=0)   #update particle positions based on available results

if (!is.null(break_flag)) break_flag=paste("nothing done; project file fulfills abortion criteria:",break_flag)


while (is.null(break_flag))
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
    update_tasklist_dds()   #update particle positions based on available results

}      
      
   
ret_val=list(value=fitness_gbest,par=X_gbest,iterations=min(iterations),break_flag=break_flag) 

return(ret_val) 
}

