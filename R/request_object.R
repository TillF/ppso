#optional functions for additional master-slave interaction

#request object from master
request_object = function(object_name, verbose_slave=FALSE)
{
  if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...requesting object '",object_name,"'"))
  mpi.send.Robj(obj=object_name, dest=0, tag=5) #tag 5 demarks request for object

  if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...waiting for object '",object_name,"'"))
  tag=0
  while(! (tag %in% c(5, 7)))                #wait till there is a message or shutdown command
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
  }
  if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  object '",object_name,"' received."))
  return(messge)
}

#push object to master
push_object = function(object_name, value, verbose_slave=FALSE)
{
  if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),": ...pushing object '",object_name,"'"))
  attr(value, "object_name")=object_name #attach name of object
  mpi.send.Robj(obj=value, dest=0, tag=6) #tag 6 demarks pushed object
  if (verbose_slave) print(paste(Sys.time(),"slave",mpi.comm.rank(),":  object '",object_name,"' sent."))
  return()
}