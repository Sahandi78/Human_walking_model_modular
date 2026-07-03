import numpy as np


class PIDController:
    """
    Joint-space PID controller with gravity
    and constraint compensation (numeric).
    """

    def __init__(self, Kp, Kd, Ki):
        self.Kp = np.diag(Kp)
        self.Kd = np.diag(Kd)
        self.Ki = np.diag(Ki)

        self.integrator_variable = np.zeros(7)

    def reset_integrator(self):
        self.integrator_variable[:] = 0.0

    def compute(
    self,
    q, dq,
    q_d, dq_d, ddq_d,
    foot_pos, dfoot_pos,
    dynamics,
    phase_name
    ):
        """
        dynamics = numeric_phases[phase_name]["dynamics"]
        """

        q = np.asarray(q).flatten()
        foot_pos = np.asarray(foot_pos).flatten()
        dq = np.asarray(dq).flatten()
        dfoot_pos = np.asarray(dfoot_pos).flatten()
        q_d = np.asarray(q_d).flatten()
        dq_d = np.asarray(dq_d).flatten()

        # --- Errors ---
        e_p = q_d - q
        e_d = dq_d - dq
        self.integrator_variable += e_p

        # --- PID ---
        tau = (
            self.Kp @ e_p +
            self.Kd @ e_d +
            self.Ki @ self.integrator_variable
        )

        # --- Gravity compensation ---
        G = dynamics["G"](np.concatenate((q, foot_pos))).flatten()
        tau += G[:7]

        # --- Double support constraint handling ---
        if phase_name.startswith("D"):
            M = dynamics["M"](np.concatenate((q, foot_pos)))
            C = dynamics["C"](np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos)))
            D = dynamics["D"](np.concatenate((q, foot_pos)))
            Gamma = dynamics["Gamma"](np.concatenate((q, foot_pos)), np.concatenate((dq, dfoot_pos))).flatten()

            q_full = np.concatenate((q, foot_pos))
            dq_full = np.concatenate((dq, dfoot_pos))

            zero = np.zeros((D.shape[0], D.shape[0]))

            Q = np.concatenate((tau, np.zeros(len(foot_pos))))
            Q -= C @ dq_full
            Q -= dynamics["G"](np.concatenate((q, foot_pos))).flatten()

            lhs = np.block([[M, D.T],
                            [D, zero]])

            rhs = np.concatenate((Q, Gamma))
            sol = np.linalg.solve(lhs, rhs)

            lambda_c = sol[len(q_full):]
            tau += D[:, :7].T @ lambda_c

        return tau.reshape(7, 1)
