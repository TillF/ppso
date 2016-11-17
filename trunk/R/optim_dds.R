optim_dds <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles= 1, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, lhc_init=FALSE, part_xchange=2,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, tryCall=FALSE, verbose=FALSE,...)
# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{

if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")
  

globvars$is_mpi = FALSE


environment(request_object)=environment() 
environment(update_tasklist_dds)=environment() #force subroutines to have this function as parent (implicitly allowing read-only access to globals)
environment(init_particles)=environment() 
environment(init_visualisation)=environment() 



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
globvars$nslaves = 1
globvars$X             =array(Inf,c(number_of_particles,number_of_parameters))  #globvars$X: position in parameter space                          
globvars$V             =array(0,c(           number_of_particles,number_of_parameters))  #globvars$V: velocity in parameter space
globvars$fitness_X     =array(Inf,number_of_particles)            
globvars$status        =array(0,number_of_particles)  #particle globvars$status: 0: to be computed; 1: finished; 2: in progress
globvars$computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if globvars$status=2)
globvars$node_id       =array(0,number_of_particles)                              #node number of worker / slave currently associated to particle
globvars$X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
globvars$fitness_lbest =array(Inf,number_of_particles)  
globvars$function_calls    =array(0,number_of_particles)  # iteration counter for each particle
                          
globvars$X_gbest     =array(Inf,number_of_parameters)            

globvars$break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the local and global fitness to the worst possible
globvars$fitness_lbest[] = Inf
globvars$fitness_gbest   = Inf

#presearch / initialisation: 
#  the particles are preferrably initialized with the data from the projectfile. If that does not exist or does not contain enough records,
#  for each uninitialized particle (uninitialized_particles) a number of prior calls (init_calls) are performed, of which the best is used
if (verbose_master) {print(paste(Sys.time(),"initializing particle positions...")); flush.console()}
  init_particles(lhc_init)  #initialize particle positions
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

    if (length(pre_run_computations) >= max_number_function_calls) stop(paste("Parameter max_number_function_calls =",max_number_function_calls,"does not suffice for initialisation. Increase it or decrease number_of_particles"))
    globvars$status[]=1; globvars$status[pre_run_computations]=0    #do computations only for the particles to be initialized, skip those that have been initialized from file
  
    globvars$computation_start[pre_run_computations]=Sys.time()
    globvars$fitness_X [pre_run_computations]=apply(as.matrix(globvars$X[pre_run_computations,],nrow=length(pre_run_computations)),1,objective_function, ...)    #execute pending pre-runs
    globvars$status    [pre_run_computations] =1      #mark as pre-runs "finished"
  
    if (!is.null(logfile))  write.table(file = logfile, cbind(format(globvars$computation_start[pre_run_computations],"%Y-%m-%d %H:%M:%S"), matrix(globvars$X[pre_run_computations,],ncol=ncol(globvars$X)), globvars$fitness_X[pre_run_computations], #write pre-runs to logfile, too
    globvars$node_id[pre_run_computations]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

   if (any(globvars$fitness_X[pre_run_computations] %in% c(NA, NaN)))
      stop("Objective function mustn't yield NA or NaN. Modify it to return very large numbers instead.")

    max_number_function_calls=max_number_function_calls-length(pre_run_computations)  #reduce number of available calls due to pre-search

    top_of_preruns=sort(globvars$fitness_X[pre_run_computations],index.return=TRUE)$ix[1:length(uninitialized_particles)]   #sort according to best fitness
    globvars$fitness_lbest[uninitialized_particles] = globvars$fitness_X [pre_run_computations[top_of_preruns]]
    globvars$X_lbest      [uninitialized_particles,] = globvars$X         [pre_run_computations[top_of_preruns],]
    calls_per_uninitialized_particle = length(pre_run_computations) %/% length(uninitialized_particles)     #distribute counting of function calls among real particles
    remaining_performed_calls = length(pre_run_computations) %% length(uninitialized_particles)
    function_calls_init = array(0,number_of_particles_org)                     #count initialisation calls extra
    function_calls_init[uninitialized_particles] = c (rep(calls_per_uninitialized_particle,   length(uninitialized_particles)-remaining_performed_calls),
                                                rep(calls_per_uninitialized_particle+1,                                 remaining_performed_calls))      
 } else
  function_calls_init = 0*globvars$function_calls

  if (verbose_master) print(paste(Sys.time()," pre-runs finished, starting actual runs..."))  
 
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

update_tasklist_dds(loop_counter=0)   #update particle positions based on available results
if (!is.null(globvars$break_flag)) globvars$break_flag=paste("nothing done; project file fulfills abortion criteria:",globvars$break_flag)


while (is.null(globvars$break_flag))
{
    globvars$computation_start[] =Sys.time()      #store time of start of this computation loop
    if (tryCall)                  #catch error message during evaluation (slower)
    {
      globvars$fitness_X=try(apply(globvars$X,1,objective_function,...),silent=TRUE)
      if (!is.numeric(globvars$fitness_X))                      #an error occured during execution
      {
        globvars$break_flag=paste("aborted: ",as.character(globvars$fitness_X))    
        next
      }        
    }
    else
      globvars$fitness_X=apply(globvars$X,1,objective_function, ...)     #no error message during evaluation (faster)
  
    globvars$status    [] =1      #mark as "finished"
    globvars$function_calls[] =globvars$function_calls[]+1        #increase iteration counter
    update_tasklist_dds()   #update particle positions based on available results

}      
      
   

ret_val=list(value=globvars$fitness_gbest,par=globvars$X_gbest, function_calls=sum(globvars$function_calls+function_calls_init),break_flag=globvars$break_flag) 

return(ret_val) 
}

