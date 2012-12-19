#create environment for global vars
#I chose global variables in order to facilitate code modularisation while keeping "passing arguments to"
# and "returning results from" functions to a minimum (speed and convenience)

globvars = new.env(hash=TRUE, parent=emptyenv())

globvars$is_mpi = (is.loaded("mpi_initialize") && mpi.comm.size()>0) #flag for distinguishing serial and parallel runs

#satisfying automatic checking
#global vars read-only (implicitly used by scoping rules, no export into globvars)
abstol = NULL
break_file = NULL
C1 = NULL
C2 = NULL
completed_particles = NULL
do_plot = NULL
do_plot_function = NULL
function_calls_init = NULL
load_projectfile = NULL
logfile = NULL
max_number_function_calls = NULL
max_number_of_iterations = NULL
max_wait_iterations = NULL
maxtries=NULL 
number_of_particles = NULL
objective_function = NULL
parameter_bounds = NULL
plot_progress = NULL
plot_window = NULL
projectfile = NULL
reltol = NULL
save_interval = NULL
tryCall = NULL
verbose_master = NULL
verbose_slave = FALSE
Vmax = NULL
w = NULL
wait_complete_iteration = NULL
wait_for_keystroke = NULL
x = NULL
y = NULL
z = NULL

#global variables with read/write access need to be imported into globvars
#globvars$break_flag
#globvars$ch
#globvars$closed_slaves
#globvars$computation_start
#globvars$evals_since_lastsave
#globvars$execution_timeout
#globvars$execution_times
#globvars$fitness_lbest
#globvars$fitness_gbest
#globvars$fitness_itbest
#globvars$fitness_X
#globvars$function_calls
#globvars$futile_iter_count
#globvars$goodness_window
#globvars$idle_slaves
#globvars$it_last_improvement
#globvars$it_last_improvement
#globvars$node_id
#globvars$slave_status
#globvars$nslaves
#globvars$relocated
#globvars$progress_window
#globvars$status
#globvars$V
#globvars$X
#globvars$X_gbest
#globvars$X_lbest


#mpi.any.source = function()
#mpi.any.tag = function()
#mpi.bcast.cmd = function()
#mpi.bcast.Robj2slave = function()
#mpi.close.Rslaves = function()
#mpi.comm.rank = function()
#mpi.comm.size = function()
#mpi.get.sourcetag = function()
#mpi.iprobe = function()
#mpi.recv.Robj = function()
#mpi.remote.exec = function()
#mpi.send.Robj = function()
#mpi.spawn.Rslaves = function()
#mpi.universe.size = function()
#slave.hostinfo = function()
#perform_task = function()
