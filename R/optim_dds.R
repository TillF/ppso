optim_dds <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles= 1, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, tryCall=FALSE)
# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{
if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")
  
eval(parse(text=paste(c("update_tasklist_dds    =",deparse(update_tasklist_dds_i )))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles         =",deparse(init_particles_i      )))))       #this creates local version of the function init_particles (see explanation there)
eval(parse(text=paste(c("init_visualisation     =",deparse(init_visualisation_i  )))))  #this creates local version of the function init_visualisation 

if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   


evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected

#initialisation
init_calls=ceiling(max(0.005*abs(max_number_function_calls),5)) #number of function calls used to initialise each particle, if no value can be loaded from a project file
number_of_particles_org=number_of_particles                 #save original number of particles
number_of_particles=max(number_of_particles,init_calls)          #increase number of particles for pre-search


nslaves = 1
X             =array(Inf,c(number_of_particles,number_of_parameters))  #X: position in parameter space                          
V             =array(0,c(           number_of_particles,number_of_parameters))  #V: velocity in parameter space
fitness_X     =array(Inf,number_of_particles)            #optimum of each particle at current iteration
status        =array(0,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
node_id       =array(0,number_of_particles)                              #node number of worker / slave
X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far
function_calls    =array(0,number_of_particles)  # iteration counter for each particle
                          
X_gbest     =array(Inf,number_of_parameters)            #global optimum

break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the local fitness to the worst possible
fitness_lbest[] = Inf
fitness_gbest = min(fitness_lbest);

#presearch / initialisation: 
#  the particles are preferrably initialized with the data from the projectfile. If that does not exist or does not contain enough records,
#  for each uninitialized particle (uninitialized_particles) a number of prior calls (init_calls) are performed, of which the best is used
  init_particles(lhc_init)  #initialize particle positions
  if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
    write.table(paste("time",paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t"),"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


if (max_number_function_calls < 0)
  {                                                         #indicator for "reset function counter" - ignore the number of function calls read from the project file
    function_calls[]=0
    max_number_function_calls=abs(max_number_function_calls)
  }
  

  status_org=status  #  store original contents

  uninitialized_particles= which(fitness_lbest[1:number_of_particles_org]==Inf)                      #"real" particles that still need to be initialized with a function value
  if (any(uninitialized_particles))
  {
    pre_run_computations = c (uninitialized_particles, seq(from=number_of_particles_org+1, length.out=max(0,init_calls - length(uninitialized_particles)))) # do preruns for uninitialized particles and the number of pending preruns

    if (length(pre_run_computations) >= max_number_function_calls) stop(paste("Parameter max_number_function_calls =",max_number_function_calls,"does not suffice for initialisation. Increase it or decrease number_of_particles"))
    status[]=1; status[pre_run_computations]=0    #do computations only for the particles to be initialized, skip those that have been initialized from file
  
    computation_start[pre_run_computations]=Sys.time()
    fitness_X [pre_run_computations]=apply(X[pre_run_computations,],1,objective_function)    #execute pending pre-runs
    status    [pre_run_computations] =1      #mark as pre-runs "finished"
  
    if (!is.null(logfile))  write.table(file = logfile, cbind(format(computation_start[pre_run_computations],"%Y-%m-%d %H:%M:%S"), matrix(X[pre_run_computations,],ncol=ncol(X)), fitness_X[pre_run_computations], #write pre-runs to logfile, too
    node_id[pre_run_computations]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

    max_number_function_calls=max_number_function_calls-length(pre_run_computations)  #reduce number of available calls due to pre-search

    top_of_preruns=sort(fitness_X[pre_run_computations],index.return=TRUE)$ix[1:length(uninitialized_particles)]   #sort according to best fitness
    fitness_lbest[uninitialized_particles] = fitness_X [pre_run_computations[top_of_preruns]]
    X_lbest      [uninitialized_particles,] = X         [pre_run_computations[top_of_preruns],]
    calls_per_uninitialized_particle = length(pre_run_computations) %/% length(uninitialized_particles)     #distribute counting of function calls among real particles
    remaining_performed_calls = length(pre_run_computations) %% length(uninitialized_particles)
    function_calls_init = function_calls                     #count initialisation calls extra
    function_calls_init[uninitialized_particles] = c (rep(calls_per_uninitialized_particle,   length(uninitialized_particles)-remaining_performed_calls),
                                                rep(calls_per_uninitialized_particle+1,                                 remaining_performed_calls))      
 } else
  function_calls_init = 0*function_calls


 
#restore array dimensions according to original number of particles
  number_of_particles=number_of_particles_org         #back to original number of particles
  X_lbest       =matrix(X_lbest      [1:number_of_particles,],ncol=number_of_parameters)        # current optimum of each particle so far
  fitness_lbest =       fitness_lbest[1:number_of_particles]                                    #best solution for each particle so far

  #restore array dimensions according to original number of particles
  X             =X_lbest  #X: position in parameter space                          
  V             =matrix(V[1:number_of_particles,],ncol=number_of_parameters)   #V: velocity in parameter space
  fitness_X     =fitness_X[1:number_of_particles]            #optimum of each particle at current iteration
#  status        =array(1,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
  computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
  node_id       =array(0,number_of_particles)                              #node number of worker / slave
  function_calls    =function_calls[1:number_of_particles]  # iteration counter for each particle
  function_calls_init = function_calls_init[1:number_of_particles]                     #count initialisation calls extra
  status=status_org[1:number_of_particles]  #  restore original contents 
  futile_iter_count = futile_iter_count[1:number_of_particles]

  fitness_gbest = min(fitness_lbest)          #update global minimum
  min_fitness_index = which(fitness_gbest==fitness_lbest)[1]
  X_gbest[] = X[min_fitness_index,]

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
    function_calls[] =function_calls[]+1        #increase iteration counter
    update_tasklist_dds()   #update particle positions based on available results

}      
      
   

ret_val=list(value=fitness_gbest,par=X_gbest,function_calls=sum(function_calls+function_calls_init),break_flag=break_flag) 

return(ret_val) 
}

