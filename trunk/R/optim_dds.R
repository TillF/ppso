optim_dds <-
function (objective_function=sample_function, number_of_parameters=2, number_of_particles=40, max_number_function_calls=500, r=0.2,  abstol=-Inf,  reltol=-Inf,  max_wait_iterations=50,
   parameter_bounds=cbind(rep(-1,number_of_parameters),rep(1,number_of_parameters)),lhc_init=FALSE,
  #runtime & display parameters
    do_plot=NULL, wait_for_keystroke=FALSE, logfile="dds.log",projectfile="dds.pro", save_interval=ceiling(number_of_particles/4),load_projectfile="try",break_file=NULL)
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

  
eval(parse(text=paste(c("update_tasklist_dds=",deparse(update_tasklist_dds_i)))))  #this creates local version of the function update_tasklist_pso (see explanation there)
eval(parse(text=paste(c("init_particles=",     deparse(init_particles_i)))))  #this creates local version of the function init_particles (see explanation there)
if ((!is.null(break_file)) && (file.exists(break_file)))      #delete break_file, if existent
  unlink(break_file)   


evals_since_lastsave=0                    #for counting function evaluations since last save of project file

# visualisation
if ((number_of_parameters!=2) | is.null(do_plot)) do_plot=FALSE           #plotting only for 2D-search
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
init_calls=ceiling(max(0.005*max_number_function_calls,5)) #number of function calls used to initialise each particle
number_of_particles_org=number_of_particles                 #save original number of particles
number_of_particles=number_of_particles*init_calls          #increase number of particles for pre-search


X             =array(0,c(number_of_particles,number_of_parameters))  #X: position in parameter space                          
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

if (!is.null(logfile) && ((load_projectfile!="loaded") || (!file.exists(logfile))))        #create logfile header, if it is not to be appended, or if it does not yet exist
  write.table(paste("time",paste(rep("parameter",number_of_parameters),seq(1,number_of_parameters),sep="_",collapse="\t"),"objective_function","worker",sep="\t") , file = logfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#presearch / initialisation
  init_particles(lhc_init)  #initialize particle positions
  number_of_particles=number_of_particles_org         #bakc to original number of particles

  #restore array dimensions according tor original number of particles
  X_lbest       =array(0.,c(number_of_particles,number_of_parameters))        # current optimum of each particle so far
  fitness_lbest =array(Inf,number_of_particles)  #best solution for each particle so far

  computation_start=Sys.time()
  fitness_X=apply(X,1,objective_function)
  status    [] =1      #mark as "finished"

  if (!is.null(logfile))  write.table(file = logfile, cbind(format(computation_start,"%Y-%m-%d %H:%M:%S"), matrix(X,ncol=ncol(X))  , fitness_X, #write pre-runs to logfile, too
  node_id), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,append=TRUE)

  max_number_function_calls=max_number_function_calls-init_calls*number_of_particles  #reduce number of available calls due to pre-search
  for (i in 1:number_of_particles)  #initialize each particle with the best of its pre-runs
  {
    offseti=(i-1)*init_calls
    min_fitness_index = which.min(fitness_X[(1:init_calls)+offseti])
    fitness_lbest[i] =fitness_X [min_fitness_index+offseti]
    X_lbest      [i,]=X         [min_fitness_index+offseti,]
  }
      
  
  #restore array dimensions according tor original number of particles
  X             =X_lbest  #X: position in parameter space                          
  V             =array(0,c(           number_of_particles,number_of_parameters))  #V: velocity in parameter space
  fitness_X     =array(Inf,number_of_particles)            #optimum of each particle at current iteration
  status        =array(0,number_of_particles)  #particle status: 0: to be computed; 1: finished; 2: in progress
  computation_start=rep(Sys.time(),number_of_particles)          #start of computation (valid only if status=2)
  node_id       =array(0,number_of_particles)                              #node number of worker / slave
  iterations    =array(0,number_of_particles)  # iteration counter for each particle
   


# actual search

fitness_itbest= fitness_gbest     #best fitness in the last it_last iterations
it_last_improvent=0               #counter for counting iterations since last improvement

while (is.null(break_flag))
{
    fitness_X=apply(X,1,objective_function)
    status    [] =1      #mark as "finished"
    iterations[] =iterations[]+1        #increase iteration counter
    update_tasklist_dds()   #update particle positions based on available results
}      
      
   
ret_val=list(value=fitness_gbest,par=X_gbest,iterations=min(iterations),break_flag=break_flag) 

return(ret_val) 
}

