from solvers.single_stance import SemiImplicitEulerSingleStance
from solvers.double_stance import SemiImplicitEulerDoubleStance

def step(t, y_current, dt, phase_name, numeric_phases, controller, q_ref, dq_ref, ddq_ref):
    """Choose correct solver based on phase type."""
    if 'D' in phase_name:
        solver = SemiImplicitEulerDoubleStance()
    else:
        solver = SemiImplicitEulerSingleStance()
    
    # Call the step method with all required arguments
    y_next, tau = solver.step(t, y_current, dt, phase_name, numeric_phases, controller, q_ref, dq_ref, ddq_ref)
    return y_next, tau