optim_pdds_robust <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=1, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, 
    tryCall=FALSE, nslaves=-1, working_dir_list=NULL, execution_timeout=NULL, verbose=FALSE)

# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{
if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")

eval(parse(text=paste(c("update_tasklist_dds    =",deparse(update_tasklist_dds_i )))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles         =",deparse(init_particles_i      )))))       #this creates local version of the function init_particles (see explanation there)
eval(parse(text=paste(c("init_visualisation     =",deparse(init_visualisation_i  )))))  #this creates local version of the function init_visualisation 
eval(parse(text=paste(c("prepare_mpi_cluster    =",deparse(prepare_mpi_cluster_i )))))  #this creates local version of the function prepare_mpi_cluster_i (see explanation there)
eval(parse(text=paste(c("check_execution_timeout=",deparse(check_execution_timeout_i)))))  #this creates local version of the function check_execution_time_i (see explanation there)
eval(parse(text=paste(c("close_mpi              =",deparse(close_mpi_i)))))  #this creates local version of the function check_execution_time_i (see explanation there)

if (!is.null(execution_timeout) && execution_timeout < 1)
{
  warning("execution_timeout must be > 1. Ignored.") 
  execution_timeout=NULL
}
  

if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   

evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected


#initialisation
init_calls=ceiling(max(0.005*abs(max_number_function_calls),5)) #number of function calls used to initialise each particle, if no value can be loaded from a project file
number_of_particles_org=number_of_particles                 #save original number of particles
number_of_particles=max(number_of_particles,init_calls)          #increase number of particles for pre-search
verbose_slave  = (verbose == TRUE) | ("slaves" %in% verbose)	#configure output level of slaves
verbose_master = (verbose == TRUE) | ("master" %in% verbose)	#configure output level of master

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

if (verbose_master) print(paste(Sys.time(),"initializing slaves..."))
if (!is.null(nslaves)) prepare_mpi_cluster(nslaves=nslaves,working_dir_list=working_dir_list,verbose_slave=verbose_slave) else nslaves=NULL             #initiate cluster, if enabled
if (verbose_master) print(paste(Sys.time(),"...slaves initialized."))

#presearch / initialisation: 
#  the particles are preferrably initialized with the data from the projectfile. If that does not exist or does not contain enough records,
#  for each uninitialized particle (uninitialized_particles) a number of prior calls (init_calls) are performed, of which the best is used
if (verbose_master) print(paste(Sys.time(),"initializing particle positions..."))
  init_particles(lhc_init)  #initialize particle positions
  if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  {
    if (!is.null(colnames(X)))
      par_names=paste(colnames(X),collapse="\t") else
      par_names=paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t") #simple numbering of parameters
    write.table(paste("time",par_names,"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
	 
	 
  if (max_number_function_calls < 0)
  {                                                         #indicator for "reset function counter" - ignore the number of function calls read from the project file
    function_calls[]=0
    max_number_function_calls=abs(max_number_function_calls)
  }

if (verbose_master)  print(paste(Sys.time(),"starting initialization runs..."))

status_org=status  #  store original contents

  uninitialized_particles= which(fitness_lbest[1:number_of_particles_org]==Inf)                      #"real" particles that still need to be initialized with a function value
  if (any(uninitialized_particles))
  {
    pre_run_computations = c (uninitialized_particles, seq(from=number_of_particles_org+1, length.out=max(0,init_calls - length(uninitialized_particles)))) # do preruns for uninitialized particles and the number of pending preruns

    if (length(pre_run_computations) >= max_number_function_calls) stop(paste("Parameter max_number_function_calls =",max_number_function_calls,"does not suffice for initialisation. Increase it or decrease number_of_particles"))
    status[]=1; status[pre_run_computations]=0    #do computations only for the particles to be initialized, skip those that have been initialized from file
  
  while ((closed_slaves < nslaves) & any(status!=1))    #do initialisation runs
  {
        tobecomputed=status==0
        while ((length(idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
        {
              current_particle=which(tobecomputed)[1]   
              slave_id=idle_slaves[1]                     #get free slave        
              if (verbose_master) print(paste(Sys.time()," ...sending task to slave",slave_id))  
              mpi.remote.exec(cmd=perform_task,params=X[current_particle,],tryCall=tryCall,slave_id=slave_id,ret=FALSE)        #submit job to slave
              if (verbose_master) print(paste(Sys.time()," ...task sent"))  
              idle_slaves=idle_slaves[-1]                         #remove this slave from list
              status            [current_particle]=2               #mark this particle as "in progress"
              node_id           [current_particle]=slave_id        #store slave_id of this task
              computation_start [current_particle]=Sys.time()      #store time of start of this computation
            tobecomputed=status==0
        } 
  
        if (!is.null(break_flag) | all(status==1)) break
        sleeptime=0
        if (verbose_master) print(paste(Sys.time()," ...wait for message from slaves..."))  
        while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #wait till there is a message
        {
#            if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
#              break_flag="user interrupt"   
            if (!is.null(execution_timeout)) sleeptime=check_execution_timeout()
            Sys.sleep(sleeptime)                                           #this prevents this loop to consume too much ressources
        }        
        if (verbose_master) print(paste(Sys.time()," ...message detected"))  
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
 
        slave_message_info <- mpi.get.sourcetag()
        slave_id <- slave_message_info[1]
        tag      <- slave_message_info[2]
        if (verbose_master) print(paste(Sys.time()," ...message received from slave",slave_id))          
      
        if (tag == 2) {      #retrieve result
          current_particle =which(node_id==slave_id & status==2)           #find which particle this result belongs to

          if (length(current_particle) > 1)                                #rr   shouldn't occur
            print(paste("strange, slave",slave_id,"returned an ambiguous result (pre-search):",slave_message))
          else    
          if (length(current_particle) ==0)                                #
          {
            if (node_interruptions[slave_id,"status"] == 0)        #rr shouldn't occur
               print(paste("strange, slave",slave_id,"returned an unrequested result (pre-search):",slave_message))          
            if (node_interruptions[slave_id,"status"] == 1)        #this is an obsolete result, from a terminated or reset slave
            {
              idle_slaves=c(idle_slaves,slave_id)     #give the slave another chance
              node_interruptions[slave_id,"status"] = 0
            }   
          } else
          {
            fitness_X [current_particle] = slave_message
            status    [current_particle] =1      #mark as "finished"
            if (!is.null(execution_timeout)) execution_times = rbind(execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),computation_start[current_particle],units="sec"))))   #monitor execution times
            if (!is.null(logfile))  write.table(file = logfile, cbind(format(computation_start[current_particle],"%Y-%m-%d %H:%M:%S"), matrix(X[current_particle,],ncol=number_of_parameters)  , fitness_X[current_particle], 
            node_id[current_particle]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)          
            idle_slaves=c(idle_slaves,slave_id)
          }
       }
        else if (tag == 3) {    # A slave has closed down.
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has gone away"))  
            closed_slaves <- closed_slaves + 1
			}
        else if (tag == 4) {    # error occured during the execution of the objective function
          if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has produced an error:",slave_message))  
		  break_flag=paste("Abort, slave",slave_id,":",slave_message)
          closed_slaves=nslaves			
        }
  }   
  
    max_number_function_calls=max_number_function_calls-length(pre_run_computations)  #reduce number of available calls due to pre-search

    top_of_preruns=sort(fitness_X[pre_run_computations],index.return=TRUE)$ix[1:length(uninitialized_particles)]   #sort according to best fitness
    fitness_lbest[uninitialized_particles] = fitness_X [pre_run_computations[top_of_preruns]]
    X_lbest      [uninitialized_particles,] = X        [pre_run_computations[top_of_preruns],]
    calls_per_uninitialized_particle = length(pre_run_computations) %/% length(uninitialized_particles)     #distribute counting of function calls among real particles
    remaining_performed_calls = length(pre_run_computations) %% length(uninitialized_particles)
    function_calls_init = function_calls                     #count initialisation calls extra
    function_calls_init[uninitialized_particles] = c (rep(calls_per_uninitialized_particle,   length(uninitialized_particles)-remaining_performed_calls),
                                                rep(calls_per_uninitialized_particle+1,                                 remaining_performed_calls))      
 } else
  function_calls_init = 0*function_calls


   if (any(fitness_X %in% c(NA, NaN)))
      stop("Objective function mustn't yield NA nro NaN. Modify it to return very large numbers instead.")

if (is.null(break_flag))
{
  if (verbose_master) print(paste(Sys.time()," pre-runs finished, starting actual runs..."))  

   #restore array dimensions according to original number of particles
  number_of_particles=number_of_particles_org         #back to original number of particles
  X_lbest       =matrix(X_lbest      [1:number_of_particles,],ncol=number_of_parameters, dimnames=dimnames(X_lbest))        # current optimum of each particle so far
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
  
  while ((closed_slaves < nslaves) )
  {
        update_tasklist_dds()   #update particle positions based on available results
        tobecomputed=status==0
  
        while ((length(idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
        {
            if (any(tobecomputed)) {
              current_particle=which.min(function_calls[tobecomputed])   #treat particles with low number of iterations first
              current_particle=which(tobecomputed)[current_particle[1]]     #choose the first entry
              slave_id=idle_slaves[1]                     #get free slave        
              if (verbose_master) print(paste(Sys.time()," ...sending task to slave",slave_id))  
              mpi.remote.exec(cmd=perform_task,params=X[current_particle,],tryCall=tryCall,slave_id=slave_id,ret=FALSE)        #submit job to slave
              if (verbose_master) print(paste(Sys.time()," ...task sent"))  
              idle_slaves=idle_slaves[-1]                         #remove this slave from list
              status            [current_particle]=2               #mark this particle as "in progress"
              node_id           [current_particle]=slave_id        #store slave_id of this task
              computation_start [current_particle]=Sys.time()      #store time of start of this computation
            }   
            if (verbose_master) print(paste(Sys.time()," ...updating task list"))  
            update_tasklist_dds()   #update particle speeds and positions based on available results
            
            tobecomputed=status==0
        } 
  
          if (!is.null(break_flag) | all(status==1)) break
  
          sleeptime=0
          if (verbose_master) print(paste(Sys.time()," ...wait for message from slaves..."))          
          while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #wait till there is a message
          {
  #            if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
  #              break_flag="user interrupt"   
              if (!is.null(execution_timeout)) sleeptime=check_execution_timeout()
              Sys.sleep(sleeptime)                                           #this prevents this loop to consume too much ressources
          }        
        if (verbose_master) print(paste(Sys.time()," ...message detected"))          
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        slave_message_info <- mpi.get.sourcetag()
        slave_id <- slave_message_info[1]
        tag      <- slave_message_info[2]
        if (verbose_master) print(paste(Sys.time()," ...message received from slave",slave_id))          
        
        if (tag == 2) {      #retrieve result
          current_particle =which(node_id==slave_id & status==2)           #find which particle this result belongs to
          if (length(current_particle) > 1)                                #rr   shouldn't occur
            print(paste("strange, slave",slave_id,"returned an ambiguous result (main search):",slave_message))
           else    
          if (length(current_particle) ==0)                                #
          {
            if (node_interruptions[slave_id,"status"] == 0)        #rr shouldn't occur
               print(paste("strange, slave",slave_id,"returned an unrequested result (main search):",slave_message))          
            if (node_interruptions[slave_id,"status"] == 1)        #this is an obsolete result, from a terminated or reset slave
            {
              idle_slaves=c(idle_slaves,slave_id)     #give the slave another chance
              node_interruptions[slave_id,"status"] = 0
            }   
          } else
          {
            fitness_X [current_particle] = slave_message
            status    [current_particle] =1      #mark as "finished"
            function_calls[current_particle] =function_calls[current_particle]+1        #increase iteration counter
            if (!is.null(execution_timeout)) execution_times = rbind(execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),computation_start[current_particle],units="sec"))))   #monitor execution times
            idle_slaves=c(idle_slaves,slave_id)
          }
        }
        else if (tag == 3) {    # A slave has closed down.
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has gone away"))  
            closed_slaves <- closed_slaves + 1
        }
        else if (tag == 4) {    # error occured during the execution of the objective function
            if (verbose_master) print(paste(Sys.time()," ...slave",slave_id,"has produced an error:",slave_message))
  		  break_flag=paste("Abort, slave",slave_id,":",slave_message)
            closed_slaves=nslaves			
        }
  }  
  if (verbose_master) print(paste(Sys.time(),"finished actual runs."))  
}

if ((closed_slaves==nslaves) && is.null(break_flag))
  break_flag = "No or delayed response from slaves" 

ret_val=list(value=fitness_gbest,par=X_gbest,function_calls=sum(function_calls+function_calls_init),break_flag=break_flag) 

if (verbose_master) print(paste(Sys.time(),"closing MPI..."))    
close_mpi()                        #diligently close Rmpi session
if (verbose_master) print(paste(Sys.time(),"...closed."))    

return(ret_val) 
}

