import numpy as np
from .base import SolverBase


class SemiImplicitEulerDoubleStance(SolverBase):
    """
    Semi-implicit Euler solver for double-stance phases
    """

    def step(self, t, y, dt, phase_name, numeric_phases, controller, q_ref, dq_ref, ddq_ref):

        # ---------------------------
        # Unpack state
        # ---------------------------
        q = y[:7]
        dq = y[7:14]
        foot_pos = y[14:18]
        dfoot_pos = y[18:22]
        lambda_prev = y[22:]

        # ---------------------------
        # Shortcuts to numeric funcs
        # ---------------------------
        numeric_phase = numeric_phases[phase_name]
        dyn = numeric_phase["dynamics"]
        con = numeric_phase["constraint"]

        M = dyn["M"]
        C = dyn["C"]
        G = dyn["G"]
        D = dyn["D"]
        Gamma = dyn["Gamma"]
        F_gen = dyn["F_generalized"]

        # ---------------------------
        # Control torques
        # ---------------------------
        tau = controller.compute(
            q=q,
            dq=dq,
            q_d=q_ref,
            dq_d=dq_ref,
            ddq_d=ddq_ref,
            foot_pos=foot_pos,
            dfoot_pos=dfoot_pos, 
            dynamics=dyn,
            phase_name=phase_name
        ).flatten()

        # ---------------------------
        # Generalized forces
        # ---------------------------
        Q = F_gen(np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos))).flatten()
        Q[:7] += tau

        # ---------------------------
        # Matrices
        # ---------------------------
        M_mat = M(np.concatenate((q, foot_pos)))
        C_mat = C(np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos)))
        G_vec = G(np.concatenate((q, foot_pos))).flatten()
        D_mat = D(np.concatenate((q, foot_pos)))
        Gamma_vec = Gamma(np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos))).flatten()

        v = np.concatenate([dq, dfoot_pos]).reshape(-1, 1)

        # ---------------------------
        # Augmented system
        # ---------------------------
        zero = np.zeros((D_mat.shape[0], D_mat.shape[0]))

        lhs = np.block([
            [M_mat, D_mat.T],
            [D_mat, zero]
        ])

        rhs = np.concatenate([
            Q - (C_mat @ v).flatten() - G_vec,
            Gamma_vec
        ])

        sol = np.linalg.solve(lhs, rhs)

        a = sol[:11]
        lambda_new = sol[11:]

        # ---------------------------
        # Semi-implicit update
        # ---------------------------
        v_new = v.flatten() + dt * a
        q_full = np.concatenate([q, foot_pos]) + dt * v_new

        # ---------------------------
        # Constraint correction (position)
        # ---------------------------
        Phi = np.vstack([
            con["Phi_1"](q_full),
            con["Phi_2"](q_full)
        ]).astype(float)
        Phi = Phi.reshape(-1, 1)

        if np.linalg.norm(Phi) > 1e-3:
            D_next = D(q_full)
            A = D_next @ D_next.T
            alpha = np.linalg.pinv(A) @ (-Phi)
            corr = D_next.T @ alpha
            q_full += corr.flatten()

        # ---------------------------
        # Velocity projection
        # ---------------------------
        M_next = M(q_full)
        D_next = D(q_full)

        Minv_Dt = np.linalg.solve(M_next, D_next.T)
        A = D_next @ Minv_Dt
        rhs_v = -(D_next @ v_new.reshape(-1, 1))

        lambda_v = np.linalg.pinv(A) @ rhs_v
        v_corr = v_new + (Minv_Dt @ lambda_v).flatten()

        # ---------------------------
        # Reassemble state
        # ---------------------------
        y_next = np.concatenate([
            q_full[:7],
            v_corr[:7],
            q_full[7:],
            v_corr[7:],
            lambda_new
        ])

        return y_next, tau
