import numpy as np

# Load symbolic phase definitions (each phase returns a dict)
from symbolic.load_phases import load_all_phases

# Substitute numerical parameters into symbolic expressions
from symbolic.substitution import substitute_phase

# Numerical values for physical / contact parameters
from symbolic.parameters import values_to_substitute

# Reference (desired) joint trajectories
from reference.desired_trajectory import DesiredTrajectory

# Utilities to Build numeric (callable) functions for each phase
from symbolic.utils.build_numeric_phase import build_numeric_phase

# PID controller imports
from controllers.pid import PIDController
from controllers.pid_instance import create_pid_controllers

# Event detection functions 
from events.gait_events import *

# Phase transition handlers
from phase_transitions.phase_layouts import PHASE_LAYOUTS
from phase_transitions.transitions import unpack_state, pack_state, transition_phase
from phase_transitions.lambda_init import initialize_lambda

# Initial condition handler
from simulation.phase_initializer import initialize_phase

# Simulation loop function
from simulation.simulation_loop import run_simulation

# Plotting results function
from plotting.tracking_plots import plot_tracking
from plotting.grf_plots import plot_grf_2x2
from plotting.torque_plots import plot_joint_torques
from plotting.animation import animate_walking

# ------------------------------------------------------------
# 1) Load ALL symbolic phase models
#    Each phase contains symbolic matrices, positions, velocities, etc.
# ------------------------------------------------------------
phases = load_all_phases()


# ------------------------------------------------------------
# 2) Substitute numerical values into each phase
#    After this step:
#      - expressions are still symbolic
#      - but all parameters (masses, lengths, etc.) are numbers
# ------------------------------------------------------------
phases_sub = {}

for name, phase in phases.items():
    phases_sub[name] = substitute_phase(phase, values_to_substitute)


# ------------------------------------------------------------
# 3) Load desired (reference) joint trajectories
#    Used later for control / tracking
# ------------------------------------------------------------
traj = DesiredTrajectory(
    q_file='reference/data/sub06_georgia/q_data.csv',
    dq_file='reference/data/sub06_georgia/dq_data.csv',
    ddq_file='reference/data/sub06_georgia/ddq_data.csv',
    dt=0.005
)


# ------------------------------------------------------------
# 4) Create numeric evaluators for ALL phases
#    This is done ONCE before the simulation loop
# ------------------------------------------------------------
numeric_phases = {}

for name, phase_sub in phases_sub.items():
    numeric_phases[name] = build_numeric_phase(phase_sub)


# ------------------------------------------------------------
# 5) Initial conditions
# ------------------------------------------------------------
start_phase = "S_L_TH"

y_current, positions, velocities = initialize_phase(
    numeric_phase=numeric_phases[start_phase],
    qd=traj.q,
    dqd=traj.dq,
    t0=0.0,
    n_lambda=4
)

# ------------------------------------------------------------
# 6) Simulation time settings
# ------------------------------------------------------------
t0 = 0.0
tf = 5.0           # total simulation time [s]
dt = 0.0001        # integration step


# ------------------------------------------------------------
# 7) Run simulation
# ------------------------------------------------------------
controller = create_pid_controllers()

results = run_simulation(
    y0=y_current,
    t0=t0,
    tf=tf,
    dt=dt,
    initial_phase=start_phase,
    numeric_phases=numeric_phases,
    controller = controller,
    q_ref_fun=traj.q,
    dq_ref_fun=traj.dq,
    ddq_ref_fun=traj.ddq
)


# ------------------------------------------------------------
# 8) Plot simulations
# ------------------------------------------------------------
plot_tracking(results)
plot_grf_2x2(results)
plot_joint_torques(results)
animate_walking(results)