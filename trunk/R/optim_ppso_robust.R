optim_ppso_robust <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40,max_number_of_iterations=5, max_number_function_calls=NULL, w=1,  C1=2, C2=2, abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   wait_complete_iteration=FALSE,parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3,lhc_init=FALSE,
  #runtime & display parameters
do_plot=NULL, wait_for_keystroke=FALSE, logfile="ppso.log",projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, tryCall=FALSE, nslaves=-1, working_dir_list=NULL, execution_timeout=NULL)
# do particle swarm optimization
{
#export r/w globals to separate environment
globvars$nslaves=nslaves
globvars$execution_timeout=execution_timeout 
 
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
#    globvars$nslaves=3                                      #number of rmpi slaves to spawn (default -1: as many as possible)
#  

if (!is.null(max_number_function_calls) && max_number_function_calls < number_of_particles)
  stop("max_number_function_calls must be at least number_of_particles.")

environment(update_tasklist_pso)=environment() #force subroutines to have this function as parent (implicitly allowing read-only access to globals)
environment(init_particles)=environment() 
environment(init_visualisation)=environment() 
environment(prepare_mpi_cluster)=environment() 
environment(check_execution_timeout)=environment() 
environment(close_mpi)=environment() 
#environment(mpi_loop)=environment() 



if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   

if (!is.null(globvars$execution_timeout) && globvars$execution_timeout < 1)
{
  warning("globvars$execution_timeout must be > 1. Ignored.") 
  globvars$execution_timeout=NULL
}

globvars$evals_since_lastsave=0                    #for counting function evaluations since last save of project file

init_visualisation()                      #prepare visualisation, if selected


#initialisation
globvars$X             =array(Inf,c(number_of_particles,number_of_parameters))  #globvars$X: position in parameter space                          
globvars$V             =array(0,c(number_of_particles,number_of_parameters))  #globvars$V: velocity in parameter space
globvars$fitness_X     =array(Inf,number_of_particles)            #optimum of each particle at current iteration
globvars$status        =array(0,number_of_particles)  #particle globvars$status: 0: to be computed; 1: finished; 2: in progress
globvars$computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if globvars$status=2)
globvars$node_id       =array(0,number_of_particles)                              #node number of worker / slave
globvars$X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
globvars$fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far
globvars$function_calls    =array(0,number_of_particles)  # globvars$function_calls counter for each particle
                          
globvars$X_gbest     =array(Inf,number_of_parameters)            #global optimum

globvars$break_flag=NULL       #flag indicating if a termination criterium has been reached


# Initialize the global and local fitness to the worst possible
 globvars$fitness_gbest = Inf;
 globvars$fitness_lbest[] = Inf

if (!is.null(globvars$nslaves)) prepare_mpi_cluster(nslaves=globvars$nslaves,working_dir_list=working_dir_list) else globvars$nslaves=NULL             #initiate cluster, if enabled

if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  {
    if (!is.null(colnames(globvars$X)))
      par_names=paste(colnames(globvars$X),collapse="\t") else
      par_names=paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t") #simple numbering of parameters
    write.table(paste("time",par_names,"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
init_particles(lhc_init)  #initialize velocities and particle positions

 
   globvars$fitness_itbest= Inf     #best fitness in the last it_last iterations
   globvars$it_last_improvement=0               #counter for counting iterations since last improvement


while ((globvars$closed_slaves < globvars$nslaves) )
{
      update_tasklist_pso()   #update particle speeds and positions based on available results
      if (!is.null(globvars$break_flag)) break
      tobecomputed=globvars$status==0
      while ((length(globvars$idle_slaves)>0) & any(tobecomputed))          #there are idle slaves available and there is work to be done
      {
          if (any(tobecomputed)) {
            current_particle=which.min(globvars$function_calls[tobecomputed])   #treat particles with low number of itereations first
            current_particle=which(tobecomputed)[current_particle[1]]     #choose the first entry
            slave_id=globvars$idle_slaves[1]                     #get free slave        
#            browser()
            mpi.remote.exec(cmd=perform_task,params=globvars$X[current_particle,],tryCall=tryCall,slave_id=slave_id,ret=FALSE)        #submit job to slave
              
            globvars$idle_slaves=globvars$idle_slaves[-1]                         #remove this slave from list
            globvars$status            [current_particle]=2               #mark this particle as "in progress"
            globvars$node_id           [current_particle]=slave_id        #store slave_id of this task
            globvars$computation_start [current_particle]=Sys.time()      #store time of start of this computation
          }   
          update_tasklist_pso()   #update particle speeds and positions based on available results
          
          tobecomputed=globvars$status==0
         flush.console()
      } 

      if (!is.null(globvars$break_flag)) break
      
      slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
      #ii: there should be some timeout here
      slave_message_info <- mpi.get.sourcetag()
      slave_id <- slave_message_info[1]
      tag      <- slave_message_info[2]
  
     
      if (tag == 2) {      #retrieve result
        current_particle =which(globvars$node_id==slave_id & globvars$status==2)           #find which particle this result belongs to
        if (length(current_particle) > 1)                                #rr   shouldn't occur
          print(paste("strange, slave",slave_id,"returned an ambiguous result (main search):",slave_message))
         else    
        if (length(current_particle) ==0)                                #
        {
          if (globvars$node_interruptions[slave_id,"status"] == 0)        #rr shouldn't occur
             print(paste("strange, slave",slave_id,"returned an unrequested result (main search):",slave_message))          
          if (globvars$node_interruptions[slave_id,"status"] == 1)        #this is an obsolete result, from a terminated or reset slave
          {
            globvars$idle_slaves=c(globvars$idle_slaves,slave_id)     #give the slave another chance
            globvars$node_interruptions[slave_id,"status"] = 0
          }   
        } else
        {
          globvars$fitness_X [current_particle] = slave_message
          globvars$status    [current_particle] =1      #mark as "finished"
          globvars$function_calls[current_particle] =globvars$function_calls[current_particle]+1        #increase iteration counter
          if (!is.null(globvars$execution_timeout)) globvars$execution_times = rbind(globvars$execution_times,data.frame(slave_id=slave_id,secs=as.numeric(difftime(Sys.time(),globvars$computation_start[current_particle],units="sec"))))   #monitor execution times
          globvars$idle_slaves=c(globvars$idle_slaves,slave_id)
        }
      }  
      else if (tag == 3) {    # A slave has closed down.
          #print(paste("slave",slave_id,"closed gracefully."))
          globvars$closed_slaves <- globvars$closed_slaves + 1
      }
      else if (tag == 4) {    # error occured during the execution of the objective function
          globvars$break_flag=paste("Abort, slave",slave_id,":",slave_message)
          globvars$closed_slaves=globvars$nslaves
      }
}      
      
if ((globvars$closed_slaves==globvars$nslaves) && is.null(globvars$break_flag))
  globvars$break_flag = "No or delayed response from slaves" 

ret_val=list(value=globvars$fitness_gbest,par=globvars$X_gbest, function_calls=sum(globvars$function_calls),break_flag=globvars$break_flag) 
  
close_mpi()                        #diligently close Rmpi session

return(ret_val) 
}

