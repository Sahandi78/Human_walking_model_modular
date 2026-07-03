from abc import ABC, abstractmethod

class SolverBase(ABC):
    """
    Base class for time integration schemes
    """

    @abstractmethod
    def step(self, t, y, dt, phase_name, numeric_phase, controller, q_ref, dq_ref, ddq_ref):
        """
        Advance the state by one time step

        Parameters
        ----------
        t : float
        y : ndarray
            Current full state
        dt : float
        phase_name : str
        phase_numeric : dict
            numeric_phases[phase_name]
        controller : PIDController

        Returns
        -------
        y_next : ndarray
        """
        pass
