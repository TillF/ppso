.onLoad <- function(libname, pkgname) {
#set MPI-flag according to startup conditions
	globvars$is_mpi=is.loaded("mpi_initialize") && mpi.comm.size()>0 #flag for distinguishing serial and parallel runs
}

