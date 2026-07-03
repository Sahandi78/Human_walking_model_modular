import numpy as np
from numpy import genfromtxt
from scipy.interpolate import interp1d


class DesiredTrajectory:
    """
    Loads desired joint trajectories and provides interpolated access
    to q(t), dq(t), ddq(t).
    """

    def __init__(self, q_file, dq_file, ddq_file, dt):
        # Load CSV data
        qd = genfromtxt(q_file, delimiter=',')
        dqd = genfromtxt(dq_file, delimiter=',')
        ddqd = genfromtxt(ddq_file, delimiter=',')

        # Transpose to shape: (n_dofs, n_time)
        self.qd = qd
        self.dqd = dqd
        self.ddqd = ddqd

        self.dt = dt
        self.n_points = self.qd.shape[1]

        self.time_span = np.arange(self.n_points) * dt

        # Interpolators
        self.q_interp = interp1d(
            self.time_span, self.qd, axis=1,
            kind='cubic', fill_value="extrapolate"
        )
        self.dq_interp = interp1d(
            self.time_span, self.dqd, axis=1,
            kind='cubic', fill_value="extrapolate"
        )
        self.ddq_interp = interp1d(
            self.time_span, self.ddqd, axis=1,
            kind='cubic', fill_value="extrapolate"
        )

    def q(self, t):
        return self.q_interp(t)

    def dq(self, t):
        return self.dq_interp(t)

    def ddq(self, t):
        return self.ddq_interp(t)
