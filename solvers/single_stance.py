import numpy as np
from .base import SolverBase


class SemiImplicitEulerSingleStance(SolverBase):

    def step(self, t, y, dt, phase_name, numeric_phases, controller, q_ref, dq_ref, ddq_ref):

        q = y[:7]
        dq = y[7:14]
        foot_pos = y[14:16]
        dfoot_pos = y[16:18]
        
        numeric_phase = numeric_phases[phase_name]
        
        dyn = numeric_phase["dynamics"]

        M = dyn["M"]
        C = dyn["C"]
        G = dyn["G"]
        F_gen = dyn["F_generalized"]

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
        
        Q = F_gen(np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos)))
        Q[:7, 0] += tau

        M_mat = M(np.concatenate((q, foot_pos)))
        C_mat = C(np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos)))
        G_vec = G(np.concatenate((q, foot_pos))).flatten()

        ddq = np.linalg.solve(
            M_mat,
            Q.flatten() - C_mat @ np.concatenate((dq, dfoot_pos)) - G_vec
        )

        dq_new = dq + dt * ddq[:7]
        dfoot_pos_new = dfoot_pos + dt * ddq[7:]
        q_new = q + dt * dq_new
        foot_pos_new = foot_pos + dt * dfoot_pos_new

        return np.concatenate([q_new, dq_new, foot_pos_new, dfoot_pos_new]), tau
