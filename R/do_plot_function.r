do_plot_function = function() #real time plotting of 2-D surface of objective function
{
    if (("base" %in% do_plot))       #do base plotting
    {
      if (exists("plot_window") && (plot_window %in% dev.list()))     #activate progress_plot window, if already open
      {
        dev.set(plot_window)
      }  else
      {
        x11()
        assign("plot_window", dev.cur(), pos=parent.frame(n=2))
      }

      res = persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",     ltheta = 120, shade = 0.75, ticktype = "detailed",      xlab = "X", ylab = "Y", zlab = "obj fun")
      points(trans3d(globvars$X[,1], globvars$X[,2], globvars$fitness_X, pmat = res), col = 2, pch =16)
      points(trans3d(globvars$X_lbest[,1], globvars$X_lbest[,2], globvars$fitness_lbest, pmat = res), col = 3, pch ="globvars$X")
      for (i in 1: min(10,number_of_particles))
         lines(rbind(trans3d(globvars$X[i,1], globvars$X[i,2], globvars$fitness_X[i], pmat = res), trans3d(globvars$X_lbest[i,1], globvars$X_lbest[i,2], globvars$fitness_lbest[i], pmat = res)) )

      if (exists("globvars$relocated") && !is.null(globvars$relocated))                   #plot particles that have been globvars$relocated during last iteration
      {
        points(trans3d(globvars$relocated[,1], globvars$relocated[,2], globvars$relocated[,3], pmat = res), col="blue", pch="O")
        assign("globvars$relocated",NULL,parent.frame(n=2))
      }

    }

    if ("rgl" %in% do_plot)        #do rgl plotting
    {
      if (!require(rgl))
      {
        warning("package rgl not found, plotting aborted")
        return()
      } 
      rgl.pop(id=hdl[c(completed_particles,completed_particles) & hdl!=0])     #remove outdated dots
      for (i in which(completed_particles))
      {
         hdl[i]                    =points3d(globvars$X[i,1],             globvars$X[i,2], globvars$fitness_X[i],     col="red")
         hdl[i+number_of_particles]=points3d(globvars$X_lbest[i,1], globvars$X_lbest[i,2], globvars$fitness_lbest[i], col="green")
      }
      if (exists("globvars$relocated") && !is.null(globvars$relocated))                   #plot particles that have been globvars$relocated during last iteration
      {
#        if (exists("hdl_r") && !is.null(hdl_r)) rgl.pop(id=hdl_r)
        hdl_r = points3d(globvars$relocated[,1], globvars$relocated[,2], globvars$relocated[,3], col="blue")
        assign("globvars$relocated",NULL,parent.frame(n=2))
        assign("hdl_r",hdl_r,parent.frame(n=2))
      }

      assign("hdl",hdl,parent.frame(n=2))
    }
}