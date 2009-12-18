#internal function: initialize particle positions and velocities

init_particles_i=function(lhc_init=FALSE)
# lhc_init: TRUE: initialise particle postions based on Latin Hypercube Sampling
#           FALSE: purely random initialisation
#note: this function reads and writes to non-local variables (i.e. varaibles declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_p*  (clumsy, but I don't know better)


{
  X=X                              #create local copies of parent variables 
  X[,]=Inf     #as a marker to denote non-initialized particles
  V=V
  X_lbest          =X_lbest
  fitness_lbest    =fitness_lbest   
  fitness_X        =fitness_X       
  status           =status          
  computation_start=computation_start
  node_id          =node_id         

  noninitialised_particles=number_of_particles      #number of particles that need to be initialized  (default:all)
  
  if ((!is.null(projectfile)) && (load_projectfile %in% c("yes","try")))
  {

    if (!file.exists(projectfile))
    {
     if (load_projectfile=="yes") stop(paste(projectfile,"could not be opened."))
    } else
    {
      proj_file_content=read.table(file = projectfile, header=TRUE,sep="\t")
      if (ncol(proj_file_content)!=3*number_of_parameters+5+1) #best_par, current_par, current_velocity)*number_of_parameters + best_objective_function+current_objective_function+status+begin_execution+node_id+function_calls
      {
        warning(paste("The number of parameters in", projectfile,"doesn't seem to match, all particles will be initialized randomly."))
#        lhc_init=FALSE
      }   else
      {
        assign("load_projectfile","loaded",                parent.frame())  #indicator that the project file has successfully been loaded
        if(!exists("number_of_particles_org")) number_of_particles_org=number_of_particles       #for DDS, the number of particles to be initialized (number_of_particles) due to the pre-run is larger than the actual number used for calculation (number_of_particles_org) 
        if (nrow(proj_file_content)>number_of_particles_org)
        {
          warning(paste(projectfile,"contains more than the specified number of",number_of_particles_org,"particles, truncated."))
          proj_file_content=proj_file_content[1:number_of_particles_org,]
        }
        noninitialised_particles=number_of_particles_org-nrow(proj_file_content)               #number of particles that need to be initialized
        if (noninitialised_particles > 0)
        {
          warning(paste(projectfile,"contains less than the specified number of",number_of_particles_org,", ",noninitialised_particles,"particle(s) will be initialized randomly."))
          noninitialised_particles=number_of_particles-nrow(proj_file_content)      #for DDS-initialisation, more particles (number_of_particles instead of number_of_particles_org) have to be initialized for the pre-run
          proj_file_content=proj_file_content[c(1:nrow(proj_file_content),rep(nrow(proj_file_content),noninitialised_particles)),names(proj_file_content)!="begin_execution"]       #just to shape the dataframe and used as marker which particles have to be initialised 
          proj_file_content[(nrow(proj_file_content)-noninitialised_particles+1):nrow(proj_file_content),]=   Inf # used as marker which particles have been initialised 
        }
  
        X_lbest           =as.matrix(proj_file_content[,1:number_of_parameters    +0])
        fitness_lbest     =as.vector(proj_file_content[,1                         +   number_of_parameters])
        X                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(1*number_of_parameters+1)])
        V                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(2*number_of_parameters+1)])
        fitness_X         =as.vector(proj_file_content[, 1                        +(3*number_of_parameters+1)])
        status            =as.vector(proj_file_content$status)
        computation_start =proj_file_content$begin_execution             #not yet used
        computation_start =strptime(computation_start,"%Y-%m-%d %H:%M:%S") #convert string to POSIX
        node_id           =as.vector(proj_file_content$node_id)
        iterations        =as.vector(proj_file_content$function_calls)
        
        node_id[status==2]=0    #any slaves marked as "in computation" in the projectfile are reset to "to be done"
        status [status==2]=0
        
        # determine the global best and its fitness from file
        min_fitness_index = which.min(fitness_lbest)
        fitness_gbest =min(fitness_lbest)          
        X_gbest[] = X_lbest[min_fitness_index[1],]
        assign("X_gbest",X_gbest,parent.frame())                         #write variable to scope of calling function
        assign("fitness_gbest",fitness_gbest,parent.frame())
        
        if (exists("Vmax") & all(V[min_fitness_index,]==0)) 
           V[min_fitness_index,]= runif(number_of_parameters,min=-0.01, max=0.01)*Vmax        #ensure that the best particle doesn't stand still
      }
    }
  }
  if (noninitialised_particles>0)          #no or not sufficient particles initialized from file -> do random initialisation
  {
    if (!("lhs" %in% installed.packages()[,"Package"]))        #check existence of lhs package
    {
      warning("Package lhs not installed, lhs_init disabled")
      lhc_init=FALSE
    }
    random_numbers=array(0,c(noninitialised_particles,number_of_parameters))
    if (!lhc_init)         #purely random initialisation
      for (i in 1 : noninitialised_particles)
        random_numbers[i,]=runif(number_of_parameters)
    else
    {
      library(lhs)
      random_numbers=randomLHS(noninitialised_particles, number_of_parameters)             #generate normalized latin hypercube realisations
    }
  
#    tobeinitialized=(number_of_particles-noninitialised_particles+1):number_of_particles      #index to particles that need to be initialized
    tobeinitialized=X[,1]==Inf      #index to particles that need to be initialized
    # Initialize the particle positions
    #X[tobeinitialized,] = parameter_bounds[,1] + (parameter_bounds[,2] - parameter_bounds[,1]) * random_numbers
    X[tobeinitialized,] = t(parameter_bounds[,1] + (parameter_bounds[,2] - parameter_bounds[,1]) * t(random_numbers)) #strangely, this transposing is necessary
    #...and their other parameters
    X_lbest       [tobeinitialized,] =  X[tobeinitialized,]
    V             [tobeinitialized,] = 0
    fitness_lbest [tobeinitialized] = Inf
    status        [tobeinitialized] = 0          
    node_id       [tobeinitialized] = 0
    iterations    [tobeinitialized] = 0
  }
  
  if (all(status==3)) status[]=0      #continue computation even if it had been finished completely before
  #"export" the variables
  assign("X",                 X,                parent.frame())
  assign("V",                 V,                parent.frame())
  assign("X_lbest",           X_lbest,          parent.frame())
  assign("fitness_lbest",     fitness_lbest,    parent.frame())
  assign("fitness_X",         fitness_X,        parent.frame())
  assign("status",            status,           parent.frame())
  assign("computation_start", computation_start,parent.frame())    #not yet used
  assign("node_id",           node_id,          parent.frame())
  assign("iterations",        iterations,       parent.frame())
}
