init_visualisation=function()
#note: this function reads and writes to non-local variables (i.e. variables declared in the calling function, usually optim_*)
#although poor style, this method was chosen to avoid passing large arrays of arguments and results, which is time-intensive
{
  # initialize visualisation
  if ((number_of_parameters!=2) | is.null(do_plot)) do_plot<<-FALSE           #plotting only for 2D-search
  if (do_plot[1]!=FALSE)
  {
    x <- seq(parameter_bounds[1,1], parameter_bounds[1,2], length= 30)
    y <- seq(parameter_bounds[2,1], parameter_bounds[2,2], length= 30)
    z=array(0,c(length(x),length(y)))
    for (i in 1: length(x))
      for (j in 1: length(y))
      z[i,j]=objective_function(c(x[i],y[j]))
    z[is.na(z)] <- 1
    
    assign("x",x,parent.frame())
    assign("y",y,parent.frame())
    assign("z",z,parent.frame())

    if ("base" %in% do_plot)  op <- par(bg = "white")       #set params for base plotting

    if ("rgl" %in% do_plot)                                 #set params for rgl plotting
    {
      if (!require(rgl))
      {
        warning("package rgl not found, rgl plotting disabled.")
        assign("do_plot",do_plot[do_plot!="rgl"],parent.frame())
      } 
      else
      {
        open3d()
        zlim <- range(y)
        zlen <- zlim[2] - zlim[1] + 1
        colorlut <- terrain.colors(zlen) # height color lookup table
        col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
        surface3d(x, y, z, color=col)
        hdl=array(0,2*number_of_particles)
        assign("hdl",hdl,parent.frame())
      }  
    }

  }

  if(is.list(plot_progress) || plot_progress==TRUE)
  {
    if (is.null(logfile) || is.null(projectfile))
    {
      warning("Output of logfile and projectfile must be enabled for plot_progress. plot_progress disabled.")
      plot_progress<<-NULL
    } else
    {
      if (is.logical(plot_progress)) plot_progress<<-NULL #remove TRUE, if  plot_progress has only been used as a binary variable
      plot_progress<<-as.list(plot_progress)
      plot_progress[["logfile"]]    <<-logfile              #prepare parameters for plot_optimization_progress
      plot_progress[["projectfile"]]<<-projectfile
    }
  } else plot_progress<<-NULL

}
