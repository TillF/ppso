optim_ppso <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40,max_number_of_iterations=5, max_number_function_calls=NULL, w=1,  C1=2, C2=2, abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   wait_complete_iteration=FALSE,parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3,lhc_init=FALSE,
  #runtime & display parameters
do_plot=NULL, wait_for_keystroke=FALSE, logfile="ppso.log",projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, nslaves=3          )
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
#    nslaves=3                                      #number of rmpi slaves to spawn (default -1: as many as possible)
# 
  
eval(parse(text=paste(c("update_tasklist_pso=",deparse(update_tasklist_pso_i)))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles=",     deparse(init_particles_i)))))  #this creates local version of the function init_particles (see explanation there)
if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   

evals_since_lastsave=0                    #for counting function evaluations since last save of project file

  prepare_mpi_cluster=function(nslaves)
  {
    if (!is.loaded("mpi_initialize")) {         
    library("Rmpi")
    }
   
    if (nslaves==-1) nslaves=mpi.universe.size() else  # Spawn as many slaves as possible
    if (nslaves>mpi.universe.size()) warning("Number of specified slaves exceeds number of available slaves.")
  
    mpi.spawn.Rslaves(nslaves=nslaves)
  
    .Last <- function(){
      if (is.loaded("mpi_initialize")){
          if (mpi.comm.size(1) > 0){
              #print("Please use mpi.close.Rslaves() to close slaves.")
              mpi.close.Rslaves()
          }
          #print("Please use mpi.quit() to quit R")
          #.Call("mpi_finalize")
      }
    }
  
    print(paste(mpi.comm.size()-1,"slaves spawned"))
    options(error=.Last)     #close rmpi on errors
  
    slavefunction <- function() {
        # Note the use of the tag for sent slave_messages:
        #     1=ready_for_task, 2=done_task, 3=exiting
        # Note the use of the tag for received slave_messages:
        #     1=task, 2=done_tasks
        junk <- 0
    
        done <- 0
        while (done != 1) {
            # Signal being ready to receive a new task
            mpi.send.Robj(junk,0,1)
    
            # Receive a task
            task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())    #this waits until a task is received
            task_info <- mpi.get.sourcetag()
            tag <- task_info[2]
    
            if (tag == 1) {
                results=task[[1]](task[[2]])
                # perform the task with the respective parameters, and create results
    
                # Send the results back as a task_done slave_message
                mpi.send.Robj(results,0,2)
                }
            else if (tag == 2) {
                # Master is saying all tasks are done.  Exit
                done <- 1
                }
            # Else ignore the slave_message or report an error
            }
    
        # Tell master that this slave is exiting.  Send master an exiting slave_message
        mpi.send.Robj(junk,0,3)
    }
    
  
    mpi.bcast.Robj2slave(objective_function)               #send function to slave
    mpi.bcast.Robj2slave(slavefunction)               #send function to slave
    
    mpi.bcast.cmd(slavefunction())        #activate all slaves by setting them to "listen" mode
    
    assign("closed_slaves",         0,parent.frame())        #set globals
    nslaves = mpi.comm.size()-1
    assign("nslaves",         nslaves,parent.frame())        #  
  }
  
  # visualisation
  if ((number_of_parameters!=2) || is.null(do_plot)) do_plot=FALSE           #plotting only for 2D-search
  if (do_plot[1]!=FALSE)
  {
    x <- seq(parameter_bounds[1,1], parameter_bounds[1,2], length= 30)
    y <- seq(parameter_bounds[2,1], parameter_bounds[2,2], length= 30)
    z=array(0,c(length(x),length(y)))
    for (i in 1: length(x))
      for (j in 1: length(y))
      z[i,j]=objective_function(c(x[i],y[j]))
    z[is.na(z)] <- 1
  
    if ("base" %in% do_plot)  op <- par(bg = "white")       #set params for base plotting
  
    if ("rgl" %in% do_plot)                                 #set params for rgl plotting  
    {
      library(rgl)
      open3d()
      zlim <- range(y)
      zlen <- zlim[2] - zlim[1] + 1
      colorlut <- terrain.colors(zlen) # height color lookup table
      col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
      surface3d(x, y, z, color=col)
      hdl=array(0,number_of_particles)
    }
  }
  
  
  #initialisation
  X             =array(0,c(number_of_particles,number_of_parameters))  #X: position in parameter space                          
  V             =array(0,c(number_of_particles,number_of_parameters))  #V: velocity in parameter space
  fitness_X     =array(Inf,number_of_particles)            #optimum of each particle at current iteration
  status        =array(0,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
  computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
  node_id       =array(0,number_of_particles)                              #node number of worker / slave
  X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
  fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far
  iterations    =array(0,number_of_particles)  # iteration counter for each particle
                            
  X_gbest     =array(Inf,number_of_parameters)            #global optimum
  
break_flag=NULL       #flag indicating if a termination criterium has been reached

# Initialize the global and local fitness to the worst possible
 fitness_gbest = Inf;
 fitness_lbest[] = Inf

init_particles(lhc_init)  #initialize velocities and particle positions

if (!is.null(logfile) & ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  write.table(paste("time",paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t"),"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

 
   fitness_itbest= Inf     #best fitness in the last it_last iterations
   it_last_improvent=0               #counter for counting iterations since last improvement

  if (!is.null(nslaves)) prepare_mpi_cluster(nslaves) else nslaves=NULL             #initiate cluster, if enabled
  
  idle_slaves=NULL  #list containing IDs of idle slaves
  junk <- 0
  
  while (closed_slaves<nslaves) 
  {
        # Receive a slave_message from a slave
#        print("receiving");    flush.console()
        slave_message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
#        print("received");    flush.console()
        #ii: there should be some timeout here
        slave_message_info <- mpi.get.sourcetag()
        slave_id <- slave_message_info[1]
        tag      <- slave_message_info[2]
    
        if (tag == 1) {         #mark slave as "ready"
            #print(paste(slave_id,"ready"))
            idle_slaves=c(idle_slaves,slave_id)
        }    
        else if (tag == 2) {      #retrieve result
          current_particle =which(node_id==slave_id & status==2)           #find which particle this result belongs to
          #ii: deal with obsolete results, deal with error message, determine average runtime
          fitness_X [current_particle] =slave_message
          status    [current_particle] =1      #mark as "finished"
          iterations[current_particle] =iterations[current_particle]+1        #increase iteration counter
          #print(paste(slave_id,"result received"))
       }
        else if (tag == 3) {    # A slave has closed down.
            #print(paste("slave",slave_id,"closed gracefully."))
            closed_slaves <- closed_slaves + 1
            if ((closed_slaves==nslaves) & (is.null(break_flag)))
                   break_flag="all slaves gone"

        }
        
        if (length(idle_slaves)>0)                 #there are idle slaves available
        {
#        print("before update");    flush.console()
           update_tasklist_pso()   #update particle speeds and positions based on available results
#           print("after update");    flush.console()
            if  (!is.null(break_flag))
            {
              slave_id=idle_slaves[length(idle_slaves)]                     #get free slave        
              idle_slaves=idle_slaves[-length(idle_slaves)]                         #remove this slave from list
#              print("before killing");    flush.console()
              mpi.send.Robj(junk, slave_id, 2)            #tell the slave to quit
#              print(paste(slave_id,"switched off")); flush.console()
              next()
            }
           tobecomputed=status==0           
           if (any(tobecomputed)) {
              current_particle=which.min(iterations[tobecomputed])   #treat particles with low number of itereations first
              current_particle=which(tobecomputed)[current_particle[1]]     #choose the first entry
              slave_id=idle_slaves[length(idle_slaves)]                     #get free slave        
              idle_slaves=idle_slaves[-length(idle_slaves)]                         #remove this slave from list
              #print("before sending");    flush.console()
              mpi.send.Robj(list(objective_function,X[current_particle,]), slave_id, 1)               # Send a task
              #print(paste("sent",db));    flush.console()
              # db=db+1
              #    if (db>199) browser()
              status            [current_particle]=2               #mark this particle as "in progress"
              node_id           [current_particle]=slave_id        #store slave_id of this task
              computation_start [current_particle]=Sys.time()      #store time of start of this computation
              #print(paste(slave_id,"command sent"))
            }    
            flush.console()
        } 
                       
        
    }
  
  ret_val=list(value=fitness_gbest,par=X_gbest,iterations=min(iterations),break_flag=break_flag) 
  
  if (!is.null(nslaves)) mpi.close.Rslaves()    #close cluster, if enabled
  return(ret_val) 
}

