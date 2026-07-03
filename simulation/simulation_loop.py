import numpy as np

from simulation.integrator import step
from simulation.grf import compute_grf
from simulation.phase_manager import update_phase
from simulation.data_recorder import DataRecorder


def run_simulation(
    y0,
    t0,
    tf,
    dt,
    initial_phase,
    numeric_phases,
    controller,
    q_ref_fun=None,
    dq_ref_fun=None,
    ddq_ref_fun=None
):
    """
    Main simulation loop.

    Parameters
    ----------
    y0 : ndarray
        Initial state vector
    t0 : float
        Initial time
    tf : float
        Final time
    dt : float
        Time step
    initial_phase : str
        Initial gait phase
    numeric_phase : dict
        Dictionary containing kinematics/velocities lambdas
    q_ref_fun, dq_ref_fun, ddq_ref_fun : callable or None
        Reference trajectory functions of time

    Returns
    -------
    results : dict
        Recorded simulation results
    """

    # --------------------------------
    # Initialization
    # --------------------------------
    t = t0
    y_current = y0.copy()
    current_phase = initial_phase

    recorder = DataRecorder(numeric_phases)
    
    torques = [0]*7

    # --------------------------------
    # Simulation loop
    # --------------------------------
    while t <= tf:

        # -----------------------------
        # Reference trajectories
        # -----------------------------
        q_ref = q_ref_fun(t) if q_ref_fun is not None else None
        dq_ref = dq_ref_fun(t) if dq_ref_fun is not None else None
        ddq_ref = ddq_ref_fun(t) if ddq_ref_fun is not None else None

        # -----------------------------
        # Compute GRFs (always 4 values)
        # -----------------------------
        F_lt, F_lh, F_rt, F_rh = compute_grf(
            t, y_current, current_phase, numeric_phases
        )


        if 'D' in current_phase:
            is_double_stance = True
        else:
            is_double_stance = False

        # -----------------------------
        # Record data BEFORE stepping
        # -----------------------------
        recorder.record(
            t=t,
            dt=dt,
            y_current=y_current,
            current_phase=current_phase,
            GRF=(F_lt, F_lh, F_rt, F_rh),
            is_double_stance=is_double_stance,
            q_ref=q_ref,
            dq_ref=dq_ref,
            ddq_ref=ddq_ref,
            torques = torques
        )

        # -----------------------------
        # Integrate dynamics
        # -----------------------------
        y_next, torques = step(t, y_current, dt, current_phase, numeric_phases, controller, q_ref, dq_ref, ddq_ref)

        # Advance time
        t_next = t + dt

        # -----------------------------
        # Update phase (events + reset)
        # -----------------------------
        y_next, phase_next = update_phase(
            t_next,
            y_next,
            current_phase,
            numeric_phases,
            F_lt, F_lh, F_rt, F_rh
        )

        # -----------------------------
        # Prepare next iteration
        # -----------------------------
        y_current = y_next
        t = t_next
        current_phase = phase_next
        torques = torques

    # --------------------------------
    # Return recorded results
    # --------------------------------
    return recorder.get_results()

