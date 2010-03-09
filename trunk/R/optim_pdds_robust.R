optim_pdds_robust <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)),lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves=-1, working_dir_list=NULL, execution_timeout=NULL)

# do Dynamically Dimensioned Search (DDS) optimization (Tolson & Shoemaker 2007)
{
if (!is.null(max_number_function_calls) && max_number_function_calls < number_of_particles)
  stop("max_number_function_calls must be at least number_of_particles.")

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
init_calls=ceiling(max(0.005*max_number_function_calls,5)) #number of function calls used to initialise each particle, if no value can be loaded from a project file
number_of_particles_org=number_of_particles                 #save original number of particles
number_of_particles=number_of_particles*init_calls          #increase number of particles for pre-search


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

if (!is.null(nslaves)) prepare_mpi_cluster(nslaves=nslaves,working_dir_list=working_dir_list) else nslaves=NULL             #initiate cluster, if enabled

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
  status_org=status  #  store original contents
  status[]=1; status[which(pre_run_computations)]=0    #do computations only for the particles to be initialized, skip those that have been initialized from file

  while ((closed_slaves < nslaves) & any(status!=1))    #do initialisation runs
  {
        tobecomputed=status==0
        while ((length(idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
        {
#            if (any(tobecomputed)) {
              current_particle=which(tobecomputed)[1]   
              slave_id=idle_slaves[1]                     #get free slave        
              mpi.remote.exec(cmd=perform_task,task=list(objective_function,X[current_particle,]),slave_id=slave_id,ret=FALSE)        #set slave to listen mode
              idle_slaves=idle_slaves[-1]                         #remove this slave from list
              status            [current_particle]=2               #mark this particle as "in progress"
              node_id           [current_particle]=slave_id        #store slave_id of this task
              computation_start [current_particle]=Sys.time()      #store time of start of this computation
#            }   
            tobecomputed=status==0
        } 
  
        if (!is.null(break_flag)) break
        sleeptime=0
        while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #wait till there is a message
        {
#        browser()
#            if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
#              break_flag="user interrupt"   
            if (!is.null(execution_timeout)) sleeptime=check_execution_timeout()
            Sys.sleep(sleeptime)                                           #this prevents this loop to consume too much ressources
        }        
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        
        slave_message_info <- mpi.get.sourcetag()
        slave_id <- slave_message_info[1]
        tag      <- slave_message_info[2]
       
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
            closed_slaves <- closed_slaves + 1
        }
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

while ((closed_slaves < nslaves) )
{
      update_tasklist_dds()   #update particle positions based on available results
      if (!is.null(break_flag)) break
      tobecomputed=status==0
      while ((length(idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
      {
          if (any(tobecomputed)) {
            current_particle=which.min(iterations[tobecomputed])   #treat particles with low number of itereations first
            current_particle=which(tobecomputed)[current_particle[1]]     #choose the first entry
            slave_id=idle_slaves[1]                     #get free slave        

            mpi.remote.exec(cmd=perform_task,task=list(objective_function,X[current_particle,]),slave_id=slave_id,ret=FALSE)        #set slave to listen mode
            idle_slaves=idle_slaves[-1]                         #remove this slave from list
            status            [current_particle]=2               #mark this particle as "in progress"
            node_id           [current_particle]=slave_id        #store slave_id of this task
            computation_start [current_particle]=Sys.time()      #store time of start of this computation
          }   
          update_tasklist_dds()   #update particle speeds and positions based on available results
          
          tobecomputed=status==0
      } 

      if (!is.null(break_flag)) break

        sleeptime=0
        while(!mpi.iprobe(mpi.any.source(),mpi.any.tag()))                #wait till there is a message
        {
#            if ((!is.null(break_file)) && (file.exists(break_file)))      #check if interrupt by user is requested
#              break_flag="user interrupt"   
            if (!is.null(execution_timeout)) sleeptime=check_execution_timeout()
            Sys.sleep(sleeptime)                                           #this prevents this loop to consume too much ressources
        }        
      slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
      slave_message_info <- mpi.get.sourcetag()
      slave_id <- slave_message_info[1]
      tag      <- slave_message_info[2]
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
          iterations[current_particle] =iterations[current_particle]+1        #increase iteration counter
          if (!is.null(execution_timeout)) execution_times = rbind(execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),computation_start[current_particle],units="sec"))))   #monitor execution times
          idle_slaves=c(idle_slaves,slave_id)
        }
      }
      else if (tag == 3) {    # A slave has closed down.
          #print(paste("slave",slave_id,"closed gracefully."))
          closed_slaves <- closed_slaves + 1
      }
      else if (tag == 4) {    # error occured during the execution of the objective function
          break_flag=paste("Abort, slave",slave_id,":",slave_message)
          closed_slaves=nslaves
      }
}      

if ((closed_slaves==nslaves) && is.null(break_flag))
  break_flag = "No or delayed response from slaves" 

ret_val=list(value=fitness_gbest,par=X_gbest,iterations=min(iterations),break_flag=break_flag) 
  
close_mpi()                        #diligently close Rmpi session

return(ret_val) 
}

