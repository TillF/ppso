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
        if (tag == 7) #this is the shutdown command
        {
            if (messge != "kill") #dunno why, but this happens!
            {
              if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  received erroneous kill message, Ignored. Message: ",messge))        
              next #wait for next message
            }
            if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  received kill message, aborting. Message: ",messge))
            assign("globvars$kill_msg", TRUE,parent.frame(n=2)) #set global flag for kill message
            return("Received kill-signal, no more waiting for requested object.")    #abort current evaluation of functions on this slave
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
     environment(get_object)=environment() 
     for (object_name in object_names) #treat multiple requests, if required
       messge[[object_name]]= get_object(object_name=object_name) #get variable or parts thereof       
     }
    
    if (length(messge)==1) messge=messge[[1]] 	  #if only one variable was requested, resolve list 
    return(messge)
  }
  
  #push object to master
  push_object = function(object_list, verbose_slave=FALSE)
  {
    if (globvars$is_mpi) #mpi-version
    {
      {
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...pushing object(s) '",paste(object_names, collapse=","),"'"))
#        attr(value, "object_name")=object_name #attach name of object
        mpi.send.Robj(obj=object_list, dest=0, tag=6) #tag 6 demarks pushed object
        if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  object '",object_name,"' sent."))
      }
    }  
    else #serial version
    {
      environment(set_object)=environment() 
      for (i in 1:length(object_list))
      {
  			object_name = names(object_list[i]) #name of object
        set_object(object_name=object_name, value=object_list[[i]]) #set variable or parts thereof
      }
    }
    return()
  }
  
set_object = function(object_name, value)     
#assign an object or parts thereof in the searchlist of environments
{
#  verbose_slave=TRUE
  if (verbose_slave)  {print(paste(Sys.time(),"setting object ",object_name));flush.console}
  
  if(!grepl("[\\[\\$]",x=object_name))     
  {
    for (env in search())
      if(exists(x=object_name, where=env))                 #try to set this object in an environment, where it already exists
      {
        assign(x=object_name, value=value, pos=env)
        break
      }
    if (!exists(x=object_name, inherits=TRUE))   #check if the object has been set somewhere
      assign(x=object_name, value=value, pos=globvars)   #assign to globvars
  } else   
  {   #only part of object requested
    if (verbose_slave) {print(paste(Sys.time(),"partial assignment"));flush.console}

    res=try(eval(parse(text=paste(object_name,"=",value)), envir = parent.frame()), silent=TRUE)
    if (class(res)=="try-error") 
      res=try(eval(parse(text=paste(object_name,"=",value)), envir = globvars), silent=TRUE)
    if (class(res)=="try-error")
    {
     if (verbose_slave) {print(paste(Sys.time(),"partial assignment unsuccessful"));flush.console}
     warning("Couldn't set value of 'object_name' on master.")
    } 
  }		  

}
  
  
get_object = function(object_name)     
#get value of specified object (or part thereof) by searching thru the environments
{  
    obj_val=NA   #return requested object, if existing, otherwise NA

    if(exists(x=object_name))
		  obj_val=get(object_name,pos=parent.frame(), inherits=FALSE) else 
	  if(exists(x=object_name, where=globvars))
		  obj_val=get(object_name,pos=globvars) else 
    if(grepl("[\\[\\$]",x=object_name))                    #only part of object requested
    {
      obj_val=try(eval(parse(text=object_name)), silent=TRUE) 
      if (class(obj_val)=="try-error") 
      {
        obj_val=try(eval(parse(text=object_name), envir = globvars), silent=TRUE) 
        if (class(obj_val)=="try-error") obj_val=NA  
      }  
    }
    return(obj_val)
}