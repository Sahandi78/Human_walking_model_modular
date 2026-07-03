import numpy as np


def initialize_lambda(t, q, dq, foot_pos, foot_vel, phase, numeric_phase):
    dyn = numeric_phase["dynamics"]

    M = dyn["M"](np.concatenate((q, foot_pos)))
    D = dyn["D"](np.concatenate((q, foot_pos)))
    C = dyn["C"](np.concatenate((q, foot_pos)), np.concatenate((dq, foot_vel)))
    G = dyn["G"](np.concatenate((q, foot_pos)))
    Gamma = dyn["Gamma"](np.concatenate((q, foot_pos)), np.concatenate((dq, foot_vel)))

    v = np.concatenate([dq, foot_vel]).reshape(-1, 1)
    Q = C @ v + G

    # Solve without explicit inverse
    Minv_Dt = np.linalg.solve(M, D.T)
    A = D @ Minv_Dt
    b = D @ np.linalg.solve(M, Q) + Gamma

    lam = np.linalg.solve(A, b)
    return lam.flatten()
