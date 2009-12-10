#test function after Ackley (1987), minimum value=-20-e
plot_optimization_progress = function  (logfile="pso.log", projectfile="pso.pro", progress_plot_filename=NULL, goodness_plot_filename=NULL, cutoff_quantile=0.95, verbose=FALSE)
{
  logfile_content    =read.table(file=  logfile,header=TRUE,sep="\t", stringsAsFactors =FALSE)
  projectfile_content=read.table(file=projectfile,header=TRUE,sep="\t", stringsAsFactors =FALSE)
  number_of_parameters=which(names(logfile_content)=="objective_function")-2

  #add positions of current minima of particles to dataset
  for (i in 1:nrow(projectfile_content))  #find at which function call number the current particle positions have been achieved
     projectfile_content$function_call_number[i]=max(-1,which(apply(apply(logfile_content[,1+(1:number_of_parameters)], 1, get("-"),t(projectfile_content[i,1:number_of_parameters]))==0,2,all)))

  if (verbose)
  {
    curbest=which.min(projectfile_content$best_objective_function)
    print(paste("current optimum found: ",projectfile_content$best_objective_function[curbest]))
    print(" at parameter set:")
    print(projectfile_content[curbest,1:number_of_parameters])
    print(paste(" found at function call ",projectfile_content$function_call_number[curbest],"from",max(projectfile_content$function_call_number),"executed calls."))
  }

  necessary_plots=number_of_parameters+2    #+2: one for objective function, one for legend
  mfrow=ceiling(sqrt(3/4*necessary_plots))    #determine number of rows and columns in plot window
  mfcol=ceiling(necessary_plots / mfrow)


#do progress plot
  if (exists("progress_window") && (progress_window %in% dev.list()))     #activate progress_plot window, if already open
  {
    dev.set(progress_window)
  }  else
  {
    x11()
    assign("progress_window", dev.cur(), pos=parent.frame())
  }

  par(mfcol=c(mfrow,mfcol))
  
  for (i in 1:number_of_parameters)
  {
    plot(logfile_content[,i+1],pch=20,xlab="function evaluations", ylab=paste("parameter",i))  
    points(projectfile_content$function_call_number,projectfile_content[,i],col="red",pch=23)
  }
  
  plot(logfile_content$objective_function,pch=20,xlab="function evaluations", ylab="objective function value",ylim=c(min(logfile_content$objective_function),quantile(logfile_content$objective_function, cutoff_quantile)),col="blue")  
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
  if (exists("goodness_window") && (goodness_window %in% dev.list()))     #activate goodness_plot window, if already open
  {
    dev.set(goodness_window)
  }  else
  {
    x11()
    assign("goodness_window", dev.cur(), pos=parent.frame())
  }
  par(mfcol=c(mfrow,mfcol))
  
  ylim=quantile(logfile_content$objective_function, c(0,cutoff_quantile))
  for (i in 1:number_of_parameters)
  {
    plot(logfile_content[,i+1],logfile_content$objective_function,pch=20,ylab="objective function value", xlab=paste("parameter",i),ylim=ylim)  
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
  for (worker in workers)
  {
    curr_worker=logfile_content$worker==worker             #get index vector to current worker
    if (length(curr_worker)<2) workers=workers[-which(workers==worker)] else   #this worker has only one call so far, omit from further treatment
    execution_time[[worker]]=diff(logfile_content$time[curr_worker])
    non_positives=execution_time[[worker]]<=0
    execution_time[[worker]][non_positives]=max(0.5,0.5*min(execution_time[[worker]][!non_positives]))    #set 0 execution times to something positive

    if (verbose)
    {
      print(paste("execution time worker",worker,"(min,median,max):"))
      cat("\t");print(min(execution_time[[worker]]))
      cat("\t");print(median(execution_time[[worker]]))
      cat("\t");print(max(execution_time[[worker]]))
    }
  }
 
  if (length(workers)>0)
    plot(range(logfile_content$time),range(unlist(execution_time)),type="n",xlab="time",ylab=paste("execution time [",attr(execution_time[[1]],"units"),"]",sep=""),log="y") #prepare plot window

  for (worker in workers)
  {
    curr_worker=which(logfile_content$worker==worker)
    points(logfile_content$time[curr_worker[-length(curr_worker)]],execution_time[[worker]],col=pal[worker],pch=(19:25)[(worker-1) %% 7 +1])
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
}