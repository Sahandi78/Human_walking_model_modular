import numpy as np


ZERO = np.zeros(3)


def compute_grf(t, y, phase_name, numeric_phases):
    """
    Compute Ground Reaction Forces (GRFs) for the current phase.

    Always returns:
        F_lt, F_lh, F_rt, F_rh

    Each force is a 2D vector [Fx, Fy].
    """
    
    # Extract q and dq
    if "S" in phase_name:
        q = np.concatenate((y[:7], y[14:16]))
        dq = np.concatenate((y[7:14], y[16:18]))
    else:
        q = np.concatenate((y[:7], y[14:18]))
        dq = np.concatenate((y[7:14], y[18:22]))
 
    # -----------------------------
    # Default: no contact anywhere
    # -----------------------------
    F_lt = ZERO.copy()
    F_lh = ZERO.copy()
    F_rt = ZERO.copy()
    F_rh = ZERO.copy()


    numeric_phase = numeric_phases[phase_name]
    contact = numeric_phase.get("contact", {})

    # Safety
    if not contact:
        return F_lt, F_lh, F_rt, F_rh

    # --------------------------------------------------
    # Single stance – LEFT
    # --------------------------------------------------
    if phase_name == "S_L_T":
        F_lt = contact["toe_force"](q, dq)

    elif phase_name == "S_L_H":
        F_lh = contact["heel_force"](q, dq)

    elif phase_name == "S_L_TH":
        F_lt = contact["toe_force"](q, dq)
        F_lh = contact["heel_force"](q, dq)

    # --------------------------------------------------
    # Double – RIGHT heel/toe entering
    # --------------------------------------------------
    elif phase_name == "D_R_H":
        F_lt = contact["left_toe_force"](q, dq)
        F_rh = contact["right_heel_force"](q, dq)

    elif phase_name == "D_R_TH":
        F_lt = contact["left_toe_force"](q, dq)
        F_rt = contact["right_toe_force"](q, dq)
        F_rh = contact["right_heel_force"](q, dq)

    # --------------------------------------------------
    # Single stance – RIGHT
    # --------------------------------------------------
    elif phase_name == "S_R_T":
        F_rt = contact["toe_force"](q, dq)

    elif phase_name == "S_R_H":
        F_rh = contact["heel_force"](q, dq)

    elif phase_name == "S_R_TH":
        F_rt = contact["toe_force"](q, dq)
        F_rh = contact["heel_force"](q, dq)

    # --------------------------------------------------
    # Double – LEFT heel/toe entering
    # --------------------------------------------------
    elif phase_name == "D_L_H":
        F_rt = contact["right_toe_force"](q, dq)
        F_lh = contact["left_heel_force"](q, dq)

    elif phase_name == "D_L_TH":
        F_rt = contact["right_toe_force"](q, dq)
        F_lt = contact["left_toe_force"](q, dq)
        F_lh = contact["left_heel_force"](q, dq)

    else:
        raise ValueError(f"Unknown phase: {phase_name}")

    return F_lt, F_lh, F_rt, F_rh
