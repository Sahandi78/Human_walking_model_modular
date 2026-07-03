import numpy as np


def initialize_state(q0, dq0, n_lambda=4):
    """
    Create initial generalized state vector y0 = [q, dq, lambdas]
    """
    lambdas0 = np.zeros(n_lambda)
    return np.concatenate((q0, dq0, lambdas0))


def initialize_kinematics(
    numeric_phase: dict,
    t0: float,
    y0: np.ndarray
):
    """
    Initialize all positions and velocities from numeric_phases[start_phase]

    Returns:
        positions: dict of lists
        velocities: dict of lists
    """

    phase = numeric_phase # A dictionary which contains lambdification functions

    positions = {}
    velocities = {}

    # -------------------------
    # Positions
    # -------------------------
    for name, (X_fun, Y_fun) in phase["positions"].items():
        positions[name] = [
            np.array([
                X_fun(t0, y0),
                Y_fun(t0, y0)
            ])
        ]

    # -------------------------
    # Velocities
    # -------------------------
    for name, V_fun in phase["velocities"].items():
        velocities[name] = [
            np.array(V_fun(t0, y0))
        ]

    return positions, velocities


def initialize_phase(
    numeric_phase: dict,
    qd,
    dqd,
    t0: float = 0.0,
    n_lambda: int = 4
):
    """
    High-level initializer used by main.py
    """

    # Generalized initial state
    y0 = initialize_state(qd(t0), dqd(t0), n_lambda)

    # Kinematics
    positions, velocities = initialize_kinematics(
        numeric_phase,
        t0,
        y0
    )

    return y0, positions, velocities
