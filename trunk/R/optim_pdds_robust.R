optim_pdds_robust <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=1, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, lhc_init=FALSE, part_xchange=2,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, 
    tryCall=FALSE, nslaves=-1, working_dir_list=NULL, execution_timeout=NULL, maxtries=10, verbose=FALSE,...)

# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{

if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")


environment(update_tasklist_dds)=environment() #force subroutines to have this function as parent (implicitly allowing read-only access to globals)
environment(init_particles)=environment() 
environment(init_visualisation)=environment() 
environment(prepare_mpi_cluster)=environment() 
environment(check_execution_timeout)=environment() 
environment(close_mpi)=environment() 
environment(mpi_loop)=environment() 

#export r/w globals to separate environment
globvars$nslaves=nslaves
globvars$execution_timeout=execution_timeout

if (!is.null(globvars$execution_timeout) && globvars$execution_timeout < 1)
{
  warning("globvars$execution_timeout must be > 1. Ignored.") 
  globvars$execution_timeout=NULL
}
  
verbose_slave  = (verbose == TRUE) | ("slaves" %in% verbose)	#configure output level of slaves
verbose_master = (verbose == TRUE) | ("master" %in% verbose)	#configure output level of master
globvars$verbose_master = verbose_master #set global value
globvars$verbose_slave  = verbose_slave  #set global value

if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   

globvars$evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected


#initialisation
init_calls=ceiling(max(0.005*abs(max_number_function_calls),5)) #number of function calls used to initialise each particle, if no value can be loaded from a project file
number_of_particles_org=number_of_particles                 #save original number of particles
#number_of_particles=max(number_of_particles,init_calls)          #increase number of particles for pre-search
number_of_particles=number_of_particles+init_calls          #increase number of particles for pre-search

#see variable explanations in globvars.R
globvars$X             =array(Inf,c(number_of_particles,number_of_parameters))  #globvars$X: position in parameter space                          
globvars$V             =array(0,c(           number_of_particles,number_of_parameters))  #globvars$V: velocity in parameter space
globvars$fitness_X     =array(Inf,number_of_particles)            
globvars$status        =array(0,number_of_particles)  #particle globvars$status: 0: to be computed; 1: finished; 2: in progress
globvars$recent_error  =array(0,number_of_particles)  # if a recent call has produced an error, the slaves's id is stored here
globvars$computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if globvars$status=2)
globvars$node_id       =array(0,number_of_particles)                              #node number of worker / slave currently associated to particle
globvars$X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
globvars$fitness_lbest =array(Inf,number_of_particles)  
globvars$function_calls    =array(0,number_of_particles)  # iteration counter for each particle
                          
globvars$X_gbest     =array(Inf,number_of_parameters)            

globvars$break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the local and global fitness to the worst possible
globvars$fitness_lbest[] = Inf
globvars$fitness_gbest = Inf

if (verbose_master) {print(paste(Sys.time(),"initializing slaves...")); flush.console()}
if (!is.null(globvars$nslaves)) prepare_mpi_cluster(nslaves=globvars$nslaves, working_dir_list=working_dir_list,verbose_slave=verbose_slave, ...) else globvars$nslaves=NULL             #initiate cluster, if enabled
if (verbose_master) {print(paste(Sys.time(),"...slaves initialized.")); flush.console()}

#presearch / initialisation: 
#  the particles are preferrably initialized with the data from the projectfile. If that does not exist or does not contain enough records,
#  for each uninitialized particle (uninitialized_particles) a number of prior calls (init_calls) are performed, of which the best is used
if (verbose_master) {print(paste(Sys.time(),"initializing particle positions...")); flush.console()}
  init_particles(lhc_init)  #initialize particle positions

if (verbose_master) {  
  if (load_projectfile=="failed")
    print(paste(Sys.time(),"...loading projectfile failed."))
  if (load_projectfile=="loaded")
    print(paste(Sys.time(),"...particles initialized from projectfile."))
  flush.console() 
}  

  if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  {
    if (!is.null(colnames(globvars$X)))
      par_names=paste(colnames(globvars$X),collapse="\t") else
      par_names=paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t") #simple numbering of parameters
    write.table(paste("time",par_names,"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
 
	 
  if (max_number_function_calls < 0)
  {                                                         #indicator for "reset function counter" - ignore the number of function calls read from the project file
    globvars$function_calls[]=0
    max_number_function_calls=abs(max_number_function_calls)
  }

if (verbose_master)  {print(paste(Sys.time(),"starting initialization runs...")); flush.console()}

globvars$status_org=globvars$status  #  store original contents

  uninitialized_particles= which(globvars$fitness_lbest[1:number_of_particles_org]==Inf)                      #"real" particles that still need to be initialized with a function value
  if (any(uninitialized_particles))
  {
    pre_run_computations = c (uninitialized_particles, seq(from=number_of_particles_org+1, length.out=max(0,init_calls - length(uninitialized_particles)))) # do preruns for uninitialized particles and the number of pending preruns

    if (length(pre_run_computations) > max_number_function_calls)
    {
     print(paste("Parameter max_number_function_calls =",max_number_function_calls,"does not suffice for initialisation. Increase it or decrease number_of_particles"))
     close_mpi()
     stop()
    } 
    globvars$status[]=1; globvars$status[pre_run_computations]=0    #do computations only for the particles to be initialized, skip those that have been initialized from file
  
    mpi_loop(init_search=TRUE, ...) #perform mpi-loop for pre-search
#    browser()
    
    max_number_function_calls=max_number_function_calls-length(pre_run_computations)  #reduce number of available calls due to pre-search

    top_of_preruns=sort(globvars$fitness_X[pre_run_computations],index.return=TRUE)$ix[1:length(uninitialized_particles)]   #sort according to best fitness
    globvars$fitness_lbest[uninitialized_particles] = globvars$fitness_X [pre_run_computations[top_of_preruns]]
    globvars$X_lbest      [uninitialized_particles,] = globvars$X        [pre_run_computations[top_of_preruns],]
    calls_per_uninitialized_particle = length(pre_run_computations) %/% length(uninitialized_particles)     #distribute counting of function calls among real particles
    remaining_performed_calls = length(pre_run_computations) %% length(uninitialized_particles)
    function_calls_init = globvars$function_calls                     #count initialisation calls extra
    function_calls_init[uninitialized_particles] = c (rep(calls_per_uninitialized_particle,   length(uninitialized_particles)-remaining_performed_calls),
                                                rep(calls_per_uninitialized_particle+1,                                 remaining_performed_calls))      
 } else
  function_calls_init = 0*globvars$function_calls


   if (any(globvars$fitness_X %in% c(NA, NaN)))
      globvars$break_flag="Objective function mustn't yield NA or NaN. Modify it to return very large numbers instead."

if (is.null(globvars$break_flag))
{
  if (verbose_master) {print(paste(Sys.time()," pre-runs finished, starting actual runs...")); flush.console()}
#  browser()
#restore array dimensions according to original number of particles
  number_of_particles=number_of_particles_org         #back to original number of particles
  globvars$X_lbest       =       globvars$X_lbest      [1:number_of_particles,,drop=FALSE]
  globvars$fitness_lbest =       globvars$fitness_lbest[1:number_of_particles]                                    #best solution for each particle so far

  #restore array dimensions according to original number of particles
  globvars$X             =globvars$X_lbest  #globvars$X: position in parameter space                          
  globvars$V             =globvars$V[1:number_of_particles,,drop=FALSE]   #globvars$V: velocity in parameter space
  globvars$fitness_X     =globvars$fitness_X[1:number_of_particles]            #optimum of each particle at current iteration
  globvars$computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if globvars$status=2)
  globvars$node_id       =array(0,number_of_particles)                              #node number of worker / slave
  globvars$function_calls    =globvars$function_calls[1:number_of_particles]  # iteration counter for each particle
  function_calls_init = function_calls_init[1:number_of_particles]                     #count initialisation calls extra
  globvars$status=globvars$status_org[1:number_of_particles]  #  restore original contents 
  globvars$status[uninitialized_particles] = 1 #formerly unitialized particles have been treated by pre-run, so set status to "finished"
  globvars$futile_iter_count = globvars$futile_iter_count[1:number_of_particles]

  globvars$fitness_gbest = min(globvars$fitness_lbest)          #update global minimum
  min_fitness_index = which(globvars$fitness_gbest==globvars$fitness_lbest)[1]
  globvars$X_gbest[] = globvars$X[min_fitness_index,]

# actual search
  if (verbose) {print("starting main search"); flush.console()}
  globvars$fitness_itbest= globvars$fitness_gbest     #best fitness in the last it_last iterations
  globvars$it_last_improvement=0               #counter for counting iterations since last improvement
  
  update_tasklist_dds()   
  if (!is.null(globvars$break_flag)) 
    globvars$break_flag=paste("nothing done; project file fulfills abortion criteria:",globvars$break_flag) else
    mpi_loop(init_search=FALSE, method="dds", ...) #perform mpi-loop for main search
 
  if (verbose_master) {print(paste(Sys.time(),"finished actual runs.")); flush.console()}
}

if ((globvars$closed_slaves==globvars$nslaves) && is.null(globvars$break_flag))
  globvars$break_flag = "No or delayed response from slaves" 

ret_val=list(value=globvars$fitness_gbest,par=globvars$X_gbest,function_calls=sum(globvars$function_calls+function_calls_init),break_flag=globvars$break_flag) 


if (verbose_master) {print(paste(Sys.time(),"closing MPI...")); flush.console()}
close_mpi()                        #diligently close Rmpi session
if (verbose_master) {print(paste(Sys.time(),"...closed.")); flush.console()}

return(ret_val) 
}

