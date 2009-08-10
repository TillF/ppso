#internal function: initialize particle positions and velocities

init_particles_i=function(lhc_init=FALSE)
# lhc_init: TRUE: initialise particle postions based on Latin Hypercube Sampling
#           FALSE: purely random initialisation
#note: this function reads and writes to non-local variables (i.e. varaibles declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
#for that purpose, this function is locally re-declared in optim_p*  (clumsy, but I don't know better)


{
  X=X                              #create local copies of parent variables 
  V=V
  X_lbest         =X_lbest
  fitness_lbest   =fitness_lbest   
  fitness_X       =fitness_X       
  status          =status          
  computation_start=computation_start
  node_id         =node_id         

  noninitialised_particles=number_of_particles      #number of particles that need to be initialized  (default:all)
  
  if ((!is.null(projectfile)) && (load_projectfile %in% c("yes","try")))
  {

    if (!file.exists(projectfile))
    {
     if (load_projectfile=="yes") stop(paste(projectfile,"could not be opened."))
    } else
    {
      proj_file_content=read.table(file = projectfile, header=TRUE,sep="\t")
      if (ncol(proj_file_content)!=3*number_of_parameters+5)
      {
        warning(paste(projectfile,"doesn't seem to match, particles will be initialized randomly."))
        lhc_init=FALSE
      }   else
      {
        assign("load_projectfile","loaded",                parent.frame())  #indicator that the project file has successfully been loaded
        if (nrow(proj_file_content)>number_of_particles)
        {
          warning(paste(projectfile,"contains more than the specified number of",number_of_particles,"particles, truncated."))
          proj_file_content=proj_file_content[1:number_of_particles,]
        }
        noninitialised_particles=number_of_particles-nrow(proj_file_content)
        if (noninitialised_particles > 0)
        {
          warning(paste(projectfile,"contains less than the specified number of",number_of_particles,"particles, missing particles will be initialized randomly."))
          proj_file_content=proj_file_content[c(1:nrow(proj_file_content),rep(nrow(proj_file_content),noninitialised_particles)),]
          proj_file_content[(nrow(proj_file_content)-noninitialised_particles+1):nrow(proj_file_content),1]=   Inf # used as marker which particles have been initialised 
          lhc_init=FALSE
        }
  
        X_lbest           =as.matrix(proj_file_content[,1:number_of_parameters    +0])
        fitness_lbest     =as.vector(proj_file_content[,1                         +   number_of_parameters])
        X                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(1*number_of_parameters+1)])
        V                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(2*number_of_parameters+1)])
        fitness_X         =as.vector(proj_file_content[, 1                        +(3*number_of_parameters+1)])
        status            =as.vector(proj_file_content[, 1                        +(3*number_of_parameters+2)])
        computation_start =proj_file_content[, 1                        +(3*number_of_parameters+3)]
        computation_start =strptime(computation_start,"%Y-%m-%d %H:%M:%S") #convert string to POSIX
        node_id           =as.vector(proj_file_content[, 1                        +(3*number_of_parameters+4)])
        
        node_id[status==2]=0    #any slaves marked as "in computation" in the projectfile are reset to "to be done"
        status [status==2]=0
        
        # determine the global best and its fitness from file
        min_fitness_index = which.min(fitness_lbest)
        fitness_gbest =min(fitness_lbest)          
        X_gbest[] = X_lbest[min_fitness_index[1],]
        assign("X_gbest",X_gbest,parent.frame())
        assign("fitness_gbest",fitness_gbest,parent.frame())
      }
    }
  }
  if (noninitialised_particles>0)          #no or not sufficient particles initialized from file -> do random initialisation
  {
    random_numbers=array(0,c(noninitialised_particles,number_of_parameters))
    if (!lhc_init)         #purely random initialisation
      for (i in 1 : noninitialised_particles)
        random_numbers[i,]=runif(number_of_parameters)
    else
    {
      library(lhs)
      random_numbers=randomLHS(noninitialised_particles, number_of_parameters)             #generate normalized latin hypercube realisations
    }
  
    tobeinitialized=(number_of_particles-noninitialised_particles+1):number_of_particles
    # Initialize the particle positions
    #X[tobeinitialized,] = parameter_bounds[,1] + (parameter_bounds[,2] - parameter_bounds[,1]) * random_numbers
    X[tobeinitialized,] = t(parameter_bounds[,1] + (parameter_bounds[,2] - parameter_bounds[,1]) * t(random_numbers)) #strangely, this transposing is necessary
    #...and their other parameters
    X_lbest       [tobeinitialized,] =  X[tobeinitialized,]
    V             [tobeinitialized,] = 0
    fitness_lbest [tobeinitialized] = Inf
    status        [tobeinitialized] = 0          
    node_id       [tobeinitialized] = 0
  }
  
  if (all(status==3)) status[]=0      #enable continuation of computation even if it had been finished properly before
  #"export" the variables
  assign("X",                 X,                parent.frame())
  assign("V",                 V,                parent.frame())
  assign("X_lbest",           X_lbest,          parent.frame())
  assign("fitness_lbest",     fitness_lbest,    parent.frame())
  assign("fitness_X",         node_id,          parent.frame())
  assign("status",            status,           parent.frame())
  assign("computation_start", computation_start,parent.frame())
  assign("node_id",           node_id,          parent.frame())
}
