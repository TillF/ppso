#plot optimization progress based on log file
#use like:
#plot_optimization_progress(logfile="ppso.log", projectfile="ppso.pro", cutoff_quantile=0.7, verbose=TRUE)

plot_optimization_progress = function  (logfile="pso.log", projectfile="pso.pro", progress_plot_filename=NULL, goodness_plot_filename=NULL, cutoff_quantile=0.95, verbose=FALSE)
{
  logfile_content    =read.table(file=  logfile,header=TRUE,sep="\t", stringsAsFactors =FALSE)
  logfile_content$objective_function=as.numeric(as.character(logfile_content$objective_function))
  projectfile_content=read.table(file=projectfile,header=TRUE,sep="\t", stringsAsFactors =FALSE)
  number_of_parameters=which(names(logfile_content)=="objective_function")-2

  #add positions of current minima of particles to dataset
  for (i in 1:nrow(projectfile_content))  #find at which function call number the current particle positions have been achieved
     projectfile_content$function_call_number[i] = 
      max(-1,which(apply(apply(logfile_content[,1+(1:number_of_parameters)], 1, get("-"),t(projectfile_content[i,1:number_of_parameters]))==0,2,all)))


  curbest=which.min(projectfile_content$best_objective_function)
  curbest_val=projectfile_content[curbest,1:number_of_parameters] #return value
  if (verbose)
  {
    print(paste("current optimum found: ",projectfile_content$best_objective_function[curbest]))
    print(" at parameter set:")
    print(projectfile_content[curbest,1:number_of_parameters])
    print(paste(" found at function call",projectfile_content$function_call_number[curbest],"of",nrow(logfile_content),"executed calls."))
  }

  necessary_plots=number_of_parameters+2    #+2: one for objective function, one for legend
  mfrow=ceiling(sqrt(3/4*necessary_plots))    #determine number of rows and columns in plot window
  mfcol=ceiling(necessary_plots / mfrow)


#do progress plot
  if (exists("progress_window", where=globvars) && (globvars$progress_window %in% dev.list()))     #activate progress_plot window, if already open
  {
    dev.set(globvars$progress_window)
  }  else
  {
    x11()
    globvars$progress_window = dev.cur()
  }

  par(mfcol=c(mfrow,mfcol))
  
  for (i in 1:number_of_parameters)
  {
    plot(logfile_content[,i+1],pch=20,xlab="function evaluations", ylab=names(logfile_content)[i+1])  
    points(projectfile_content$function_call_number,projectfile_content[,i],col="red",pch=23)
  }
  
  ylim=quantile(logfile_content$objective_function[is.finite(logfile_content$objective_function)],
                c(0,cutoff_quantile))
  plot(logfile_content$objective_function,pch=20,xlab="function evaluations", ylab="objective function value",ylim = ylim,col="blue")  
  points(projectfile_content$function_call_number,projectfile_content$best_objective_function,col="red",pch=23)
  plot.new()
  legend("left", "particles' best", pch=23, col="red")    #plot legend only once
  if (!is.null(progress_plot_filename))
  {
    tt=regexpr("[^\\.]*$",progress_plot_filename)
    extension=substring(progress_plot_filename,tt)
    if (!(extension %in% c("wmf", "emf", "png", "jpg", "jpeg", "bmp", "tif", "tiff", "ps", "eps", "pdf")))
    {
      extension="ps"
      progress_plot_filename=paste(progress_plot_filename,".ps",sep="")
    }
    savePlot(filename = progress_plot_filename, type = extension, device = dev.cur())
  }

#do goodness plot
  if (exists("goodness_window", where=globvars) && (globvars$goodness_window %in% dev.list()))     #activate goodness_plot window, if already open
  {
    dev.set(globvars$goodness_window)
  }  else
  {
    x11()
    globvars$goodness_window=dev.cur()
  }
  par(mfcol=c(mfrow,mfcol))
  
  for (i in 1:number_of_parameters)
  {
    plot(logfile_content[,i+1],logfile_content$objective_function,pch=20,ylab="objective function value", xlab=names(logfile_content)[i+1], ylim=ylim)  
    points(projectfile_content[,i],projectfile_content$best_objective_function,col="red",pch=23)
  }
  workers=sort(unique(logfile_content$worker))
  if (length(workers)>=2)
    pal=palette(topo.colors(length(workers))) else  #set palette
  {
    pal="blue"  #single worker application, treat differently
    workers=1                 
    logfile_content$worker=1
  }
  
  # plot workload distribution
  hist(logfile_content$worker,plot=TRUE,xlab="worker ID",ylab="#function calls",breaks=seq(min(logfile_content$worker)-0.5,max(logfile_content$worker)+0.5,by=1),col=pal,main="")
#  #plot worker history
  logfile_content$time=as.POSIXct(strptime(logfile_content$time,format="%Y-%m-%d %H:%M:%S"))
#  plot(logfile_content$time,logfile_content$worker,main="worker history", xlab="time",ylab="worker ID")
#  
# plot execution time history
  execution_time=list()
  tunits=NULL
  omitted_workers=NULL
  mintime=Inf
  for (i in 1:length(workers))
  {
    worker=workers[i]
    curr_worker=logfile_content$worker==worker             #get index vector to current worker
    ncalls = sum(curr_worker)
    if (ncalls < 2) 
    {
      omitted_workers=c(omitted_workers, worker)    #this worker has only one call so far, omit from further treatment
      execution_time[[i]] = NA
      next
    } else
    execution_time[[i]]=diff(logfile_content$time[curr_worker])

    
    positive_times=execution_time[[i]]>0
    mintime=min(c(mintime, execution_time[[i]][positive_times]))
   
    tunits=c(tunits,attr(execution_time[[i]],"units")) #collect units (may differ between slaves)
  }

  #ensure same units in time
  most_used_unit = names(sort(table(tunits), decreasing = TRUE)[1]) #find most used units among execution times
  for (i in 1:length(workers))
  { 
    units(execution_time[[i]])=most_used_unit
    
    if (verbose)
    {
      print(paste("execution time worker",worker,"(min,median,max):"))
      cat("\t");print(min   (execution_time[[i]]))
      cat("\t");print(median(execution_time[[i]]))     
      cat("\t");print(max   (execution_time[[i]]))
    }
    
    positive_times=execution_time[[i]]>0
    execution_time[[i]][!positive_times]=mintime    #set execution times that are zero to something positive to allow log plot
    
  }  
  if (mintime == 0)
    logplot = "" else
    logplot = "y"
  

  if (length(workers)- length(omitted_workers) > 0)
  {
    plot(range(logfile_content$time),range(unlist(execution_time), na.rm=TRUE),type="n",xlab="time",ylab=paste("execution time [",attr(execution_time[[1]],"units"),"]",sep=""),log=logplot) #prepare plot window

    for (i in 1:length(workers))
    {
      worker=workers[i]
      if (worker %in% omitted_workers) next
      curr_worker=which(logfile_content$worker==worker)
      points(logfile_content$time[curr_worker[-length(curr_worker)]],execution_time[[i]],col=pal[worker],pch=".") #pch=(19:25)[(worker-1) %% 7 +1]
    }
  }  

  if (!is.null(goodness_plot_filename))
  {
    tt=regexpr("[^\\.]*$",goodness_plot_filename)
    extension=substring(goodness_plot_filename,tt)
    if (!(extension %in% c("wmf", "emf", "png", "jpg", "jpeg", "bmp", "tif", "tiff", "ps", "eps", "pdf")))
    {
      extension="ps"
      goodness_plot_filename=paste(goodness_plot_filename,".ps",sep="")
    }
    savePlot(filename = goodness_plot_filename, type = extension, device = dev.cur())
  }
  
  return(curbest_val)
}