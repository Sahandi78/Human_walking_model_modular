import numpy as np
from .phase_layouts import PHASE_LAYOUTS
from .lambda_init import initialize_lambda

def unpack_state(y, phase_name):
    layout = PHASE_LAYOUTS[phase_name]
    
    q = y[:7]
    dq = y[7:14]

    idx = 14
    n_feet = layout["n_feet"]

    foot_pos = y[idx : idx + 2*n_feet]
    idx += 2*n_feet

    foot_vel = y[idx : idx + 2*n_feet]
    idx += 2*n_feet

    if layout["has_lambda"]:
        lam = y[idx:]
    else:
        lam = None

    return q.copy(), dq.copy(), foot_pos.copy(), foot_vel.copy(), lam




def pack_state(q, dq, foot_pos, foot_vel, phase_name, lam=None):
    layout = PHASE_LAYOUTS[phase_name]

    y = [q, dq, foot_pos, foot_vel]

    if layout["has_lambda"]:
        if lam is None:
            raise ValueError("Lambda required for double-stance phase")
        y.append(lam)

    return np.concatenate(y)




def transition_phase(
    t,
    y,
    from_phase,
    to_phase,
    numeric_phases,
    new_contacts
):
    """
    new_contacts: dict like
    {
        "R_heel": (X, Y, dX, dY),
        "L_toe": (X, Y, dX, dY)
    }
    """

    q, dq, foot_pos_old, foot_vel_old, _ = unpack_state(y, from_phase)

    # Assemble new foot vectors in correct order
    foot_pos = []
    foot_vel = []

    for label in PHASE_LAYOUTS[to_phase]["foot_labels"]:
        X, Y, dX, dY = new_contacts[label]
        foot_pos.extend([X, Y])
        foot_vel.extend([dX, dY])

    foot_pos = np.asarray(foot_pos)
    foot_vel = np.asarray(foot_vel)

    # Initialize lambda if needed
    if PHASE_LAYOUTS[to_phase]["has_lambda"]:
        lam = initialize_lambda(
            t, q, dq, foot_pos, foot_vel, to_phase, numeric_phases[to_phase]
        )
    else:
        lam = None

    return pack_state(q, dq, foot_pos, foot_vel, to_phase, lam)
