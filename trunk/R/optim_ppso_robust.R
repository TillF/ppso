optim_ppso_robust <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40,max_number_of_iterations=5, max_number_function_calls=500, w=1,  C1=2, C2=2, abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   wait_complete_iteration=FALSE,parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)), initial_estimates=NULL, Vmax=(parameter_bounds[,2]-parameter_bounds[,1])/3,lhc_init=FALSE,
  #runtime & display parameters
do_plot=NULL, wait_for_keystroke=FALSE, logfile="ppso.log",projectfile="ppso.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL, plot_progress=FALSE, 
tryCall=FALSE, nslaves=-1, working_dir_list=NULL, execution_timeout=NULL, maxtries=10, verbose=FALSE, ...)
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

if (!is.null(max_number_function_calls) && abs(max_number_function_calls) < number_of_particles)
  stop("abs(max_number_function_calls) must be at least number_of_particles.")

environment(update_tasklist_pso)=environment() #force subroutines to have this function as parent (implicitly allowing read-only access to globals)
environment(init_particles)=environment() 
environment(init_visualisation)=environment() 
environment(prepare_mpi_cluster)=environment() 
environment(check_execution_timeout)=environment() 
environment(close_mpi)=environment() 
environment(mpi_loop)=environment() 

verbose_slave  = (verbose == TRUE) | ("slaves" %in% verbose)	#configure output level of slaves
verbose_master = (verbose == TRUE) | ("master" %in% verbose)	#configure output level of master
globvars$verbose_master = verbose_master #set global value
globvars$verbose_slave  = verbose_slave  #set global value


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
globvars$recent_error  =array(0,number_of_particles)  # if a recent call has produced an error, the slaves's id is stored here
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

if (verbose_master) {print(paste(Sys.time(),"initializing slaves...")); flush.console()}
if (!is.null(globvars$nslaves)) prepare_mpi_cluster(nslaves=globvars$nslaves, working_dir_list=working_dir_list,verbose_slave=verbose_slave, ...) else globvars$nslaves=NULL             #initiate cluster, if enabled
if (verbose_master) {print(paste(Sys.time(),"...slaves initialized.")); flush.console()}


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
   if (verbose_master) {print(paste(Sys.time(),"...prepared log file (header written)")); flush.console()}
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

   globvars$fitness_itbest= Inf     #best fitness in the last it_last iterations
   globvars$it_last_improvement=0               #counter for counting iterations since last improvement

 update_tasklist_pso()   
  if (!is.null(globvars$break_flag)) 
    globvars$break_flag=paste("nothing done; project file fulfills abortion criteria:",globvars$break_flag) else
    mpi_loop(init_search=FALSE, method="pso", ...) #perform mpi-loop 
 
  if (verbose_master) print(paste(Sys.time(),"finished actual runs."))  


if ((globvars$closed_slaves==globvars$nslaves) && is.null(globvars$break_flag))
  globvars$break_flag = "No or delayed response from slaves" 

ret_val=list(value=globvars$fitness_gbest,par=globvars$X_gbest,function_calls=sum(globvars$function_calls),break_flag=globvars$break_flag) 


if (verbose_master) {print(paste(Sys.time(),"closing MPI...")); flush.console()}
close_mpi()                        #diligently close Rmpi session
if (verbose_master) {print(paste(Sys.time(),"...closed.")); flush.console()}

return(ret_val) 
}

