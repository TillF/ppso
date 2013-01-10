#optional functions for additional master-slave interaction

#mpi-versions
  #request object from master
  request_object = function(object_names, verbose_slave=FALSE)
  {
    if (globvars$is_mpi) #mpi-version
    {
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...requesting object(s) '",paste(object_names, collapse=",") ,"'"))
      mpi.send.Robj(obj=object_names, dest=0, tag=5) #tag 5 demarks request for object
    
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...waiting for object(s) '",paste(object_names, collapse=","),"'"))
      tag=0
      while(! (tag %in% c(5, 7)))                #wait till there is a reply message or shutdown command
      {
        messge = mpi.recv.Robj(source=0, tag=mpi.any.tag())
        messge_info = mpi.get.sourcetag()
        tag      = messge_info[2]
        if (tag == 7)
        {
            if (messge != "kill") #dunno why, but this happens!
            {
              if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  received erroneous kill message, Ignored. Message: ",messge))        
              next
            }
            if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  received kill message, aborting. Message: ",messge))
            stop("slave stopped.")    #abort current evaluation of functions on this slave
            return("xstopped")
        }
        if (is.function(messge) || length(messge)!=length(object_names)) 
        {
          tag=2
          next #strangely, this occurs, i.e. sometime wrongly tagged messages are received
        }  
          
      }
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  object '",object_names,"' received."))
    }  else   #serial version
    {         
     if (verbose_slave) print(paste(Sys.time(),": ...requesting object '",paste(object_names, collapse=","),"'"))

     messge=list()
     for (object_name in object_names) #treat multipiple requests, if required
     {
       if(exists(x=object_name))
        messge[[object_name]]=get(object_name,pos=parent.frame(n=2), inherits=FALSE) else 
       if(exists(x=object_name, where=globvars))
        messge[[object_name]]=get(object_name,pos=globvars) else   #return requested object, if existing, otherwise NA
    	  messge[[object_name]]=NA
     }
    }
    
    if (length(messge)==1) messge=messge[[1]] 	  #if only one variable was requested, resolve list 
    return(messge)
  }
  
  #push object to master
  push_object = function(object_name, value, verbose_slave=FALSE)
  {
    if (globvars$is_mpi) #mpi-version
    {
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...pushing object '",object_name,"'"))
      attr(value, "object_name")=object_name #attach name of object
      mpi.send.Robj(obj=value, dest=0, tag=6) #tag 6 demarks pushed object
      if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  object '",object_name,"' sent."))
    } else #serial version
    {
      for (env in search())
        if(exists(x=object_name, where=env))                 #try to set this object in an environment, where it already exists
        {
          assign(x=object_name, value=value, pos=env)
          break
        }  
      if (!exists(x=object_name, inherits=TRUE))   #check if the object has not been found somewhere
        assign(x=object_name, value=value, pos=globvars)    
    }
    return()
  }