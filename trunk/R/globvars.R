#create environment for global vars
#I chose global variables in order to facilitate code modularisation while keeping "passing arguments to"
# and "returning results from" functions to a minimum (speed and convenience)

#globvars = new.env(hash=TRUE, parent=emptyenv())
globvars = new.env(hash=TRUE)

#satisfying automatic checking
#global vars read-only (implicitly used by scoping rules, no export into globvars)
abstol = NULL
break_file = NULL
C1 = NULL
C2 = NULL
completed_particles = NULL
do_plot = NULL
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
#globvars$break_flag  #verbose indicator for interruption status
#globvars$ch        #character for reading inout in interactive mode 
#globvars$ch_mpi    #character for reading inout in interactive mode (MPI-Loop)
#globvars$closed_slaves
#globvars$computation_start
#globvars$evals_since_lastsave
#globvars$execution_timeout
#globvars$execution_times
#globvars$fitness_lbest       #best solution for each particle so far
#globvars$fitness_gbest       #current global optimum
#globvars$fitness_itbest
#globvars$fitness_X           #fitness of each particle at current iteration
#globvars$function_calls
#globvars$futile_iter_count
#globvars$goodness_window
#globvars$idle_slaves
#globvars$it_last_improvement
#globvars$it_last_improvement
#globvars$mpi_mode              #"bcast": slave are idle, tasks are broadcasted to all and attended only by one;"loop": slaves run in loop waiting for message, tasks are send to specific slaves only
#globvars$node_id               #[1:number_of_particles] #node number of worker / slave currently associated to particle
#globvars$slave_status	[1:nslaves]
#               "counter"           #counts total number of timeouts produced by slave
#               "timeouts_in_row"   #counts successive timeouts
#               "task_particle"     # holds particle-ID that is currently treated
#globvars$nslaves
#globvars$relocated
#globvars$progress_window
#globvars$status                    # [1:nparticles]: denotes open task (0), finished (1) or computing (2)
#globvars$recent_error                    # if a recent call has produced an error, the slaves's id is stored here

#globvars$V
#globvars$X                 #current positions of particles (number_of_particles X number_of_parameters)
#globvars$X_gbest           #position of global best solution so far [1:number_of_parameters]
#globvars$X_lbest           #position of current optimum of each particle so far  (number_of_particles X number_of_parameters)


#do_plot_function = function()
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
