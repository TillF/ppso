#internal function: initialize particle positions and velocities

init_particles=function(lhc_init=FALSE)
# lhc_init: TRUE: initialise particle postions based on Latin Hypercube Sampling
#           FALSE: purely random initialisation
#note: this function reads and writes to non-local variables (i.e. varaibles declared in the calling function, usually optim_p*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive

{
  globvars$X[,]=Inf     #as a marker to denote non-initialized particles

  param_names=rownames(parameter_bounds)     #try to retrieve parameter names
  if (length(param_names)==0)
    param_names=rownames(initial_estimates)
  
  if (length(param_names)!=0)
  {
    colnames(globvars$X)=param_names
    colnames(globvars$X_lbest)=param_names
    names(globvars$X_gbest)=param_names
  }
  noninitialised_particles=number_of_particles      #number of particles that need to be initialized  (default:all)
  if(!exists("number_of_particles_org",parent.frame(), inherits=FALSE)) number_of_particles_org=number_of_particles       #for DDS, the number of particles to be initialized (number_of_particles) due to the pre-run is larger than the actual number used for calculation (number_of_particles_org) 
  
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
          proj_file_content=proj_file_content[c(1:nrow(proj_file_content),rep(nrow(proj_file_content),noninitialised_particles)), names(proj_file_content) != "begin_execution"]       #just to shape the dataframe and used as marker which particles have to be initialised 
          proj_file_content[(nrow(proj_file_content)-noninitialised_particles+1):nrow(proj_file_content),]=   Inf # used as marker which particles have been initialised 
          proj_file_content$begin_execution=ISOdate(1970,1,1)+NA #discard begin_execution, as it is no longer meaningful
        }
  
        globvars$X_lbest           =as.matrix(proj_file_content[,1:number_of_parameters    +0])
        globvars$fitness_lbest     =as.vector(proj_file_content[,1                         +   number_of_parameters])
        globvars$X                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(1*number_of_parameters+1)])
        globvars$V                 =as.matrix(proj_file_content[,(1:number_of_parameters)  +(2*number_of_parameters+1)])
        globvars$fitness_X         =as.vector(proj_file_content[, 1                        +(3*number_of_parameters+1)])
#        globvars$status            =as.vector(globvars$status)
        globvars$computation_start =proj_file_content$begin_execution             
        globvars$computation_start =strptime(globvars$computation_start,"%Y-%m-%d %H:%M:%S") #convert string to POSIX
#        globvars$node_id           =as.vector(globvars$node_id)
#        globvars$function_calls    =as.vector(globvars$function_calls)
        globvars$function_calls    =as.vector(proj_file_content$function_calls)
        
        globvars$node_id[globvars$status==2]=0    #any slaves marked as "in computation" in the projectfile are reset to "to be done"
        globvars$status [globvars$status==2]=0
        
        min_fitness_index = which.min(globvars$fitness_lbest)
        if (exists("Vmax") & all(globvars$V[min_fitness_index,]==0)) 
           globvars$V[min_fitness_index,]= runif(number_of_parameters,min=-0.01, max=0.01)*Vmax        #ensure that the best particle doesn't stand still
      }
    }
  }

  if (!is.null(initial_estimates))
  {
    initial_estimates = as.matrix(initial_estimates)       #reshape any vector as matrix
    if (nrow(initial_estimates)!=number_of_parameters)
    {
      warning("initial_estimates must contain <number_of_parameters> rows, ignored.") 
      initial_estimates =NULL
    }  else
    {
      out_of_bounds=NULL    #initial estimates that are out of bounds
      for (i in 1:ncol(initial_estimates))
        if (any ((initial_estimates[,i] < parameter_bounds[,1]) | (initial_estimates[,i] > parameter_bounds[,2])))
         out_of_bounds=c(out_of_bounds,i)
      if (any(out_of_bounds))
      {
        warning(paste("initial estimates",paste(out_of_bounds,collapse=", ")," out of bounds, ignored"))
        initial_estimates = initial_estimates[,- out_of_bounds]     #discard invalid initial estimates
      }
      if (ncol(initial_estimates) > noninitialised_particles)
      {
        warning(paste ("sufficient initial estimates loaded from project file, ", ncol(initial_estimates) - noninitialised_particles, "columns of argument <initial_estimates> ignored"))
        initial_estimates = initial_estimates[,1:noninitialised_particles]     #discard invalid initial estimates
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
    tobeinitialized = globvars$X[,1]==Inf      #index to particles that need to be initialized
    # Initialize the particle positions
    globvars$X[tobeinitialized,] = t(parameter_bounds[,1] + (parameter_bounds[,2] - parameter_bounds[,1]) * t(random_numbers)) #strangely, this transposing is necessary
    if (!is.null(initial_estimates) && (ncol(initial_estimates)>0))         #if any initial estimates have been specified as an argument, use these
      globvars$X[which(tobeinitialized)[1:ncol(initial_estimates)],] = t(initial_estimates)

    #...and their other parameters
    globvars$X_lbest       [tobeinitialized,] =  globvars$X[tobeinitialized,]
    globvars$V             [tobeinitialized,] = 0
    globvars$fitness_lbest [tobeinitialized] = Inf
    globvars$status        [tobeinitialized] = 0          
    globvars$node_id       [tobeinitialized] = 0
    globvars$function_calls    [tobeinitialized] = 0
  }
  
  # determine the global best and its fitness from all available data
  min_fitness_index = which.min(globvars$fitness_lbest)
  globvars$fitness_gbest =min(globvars$fitness_lbest)          
  globvars$X_gbest[] = globvars$X_lbest[min_fitness_index[1],]

  if (length(param_names)!=0)
  {
    colnames(globvars$X)=param_names
    colnames(globvars$X_lbest)=param_names
    names(globvars$X_gbest)=param_names
  }
  
    
  if (all(globvars$status==3)) globvars$status[]=0      #continue computation even if it had been finished completely before
  
  globvars$futile_iter_count=array(0,number_of_particles_org)
  
  globvars$slave_status=array(0,c(globvars$nslaves,2),dimnames=list(NULL,c("counter","timeouts_in_row")))         #array for recording slave fitness
  if (!exists("globvars$execution_timeout") || !is.null(globvars$execution_timeout))                          #monitor execution times
    globvars$execution_times=data.frame(slave_id=NULL,secs=NULL)
    
}
