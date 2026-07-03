import sympy as sp
import numpy as np



#-------- Position lambdification function 
def lambdify_position(expr, phase):
    """
    Lambdify a scalar position expression for a given phase.
    """
    q_syms = phase["q"]
    lam = sp.lambdify(q_syms, expr, modules="numpy")

    if len(q_syms) == 9:
        def f(t, y):
            q = np.concatenate((y[:7], y[14:16]))
            return lam(*q)
    else:
        def f(t, y):
            q = np.concatenate((y[:7], y[14:18]))
            return lam(*q)

    return f



#-------- Velocity lambdification function 
def lambdify_velocity(expr, phase):
    """
    Lambdify a velocity expression for a given phase.
    """
    q_syms = phase["q"]
    dq_syms = phase["dq"]

    lam = sp.lambdify(q_syms + dq_syms, expr, modules="numpy")

    if len(q_syms) == 9:
        def f(t, y):
            q = np.concatenate((y[:7], y[14:16]))
            dq = np.concatenate((y[7:14], y[16:18]))
            return lam(*(np.concatenate((q, dq))))
    else:
        def f(t, y):
            q = np.concatenate((y[:7], y[14:18]))
            dq = np.concatenate((y[7:14], y[18:22]))
            return lam(*(np.concatenate((q, dq))))

    return f
    


#-------- Dynamics lambdification function 
def lambdify_matrix(expr, phase_sub, use_dq=False):
    """
    Lambdify a matrix-valued SymPy expression for a given phase.
    Parameters
    ----------
    expr : sympy.Matrix
        Symbolic matrix after constant substitution
    phase_sub : dict
        Phase dictionary containing 'q' and 'dq'
    use_dq : bool
        Whether the matrix depends on dq
    Returns
    -------
    function
        f(q) or f(q, dq) → numpy array
    """
    q_syms = phase_sub["q"]
    dq_syms = phase_sub["dq"]
    
    if use_dq:
        args = q_syms + dq_syms
        lam = sp.lambdify(args, expr, modules="numpy")
        
        def f(q, dq):
            return lam(*(np.concatenate((q, dq))))
        return f
    else:
        args = q_syms
        lam = sp.lambdify(args, expr, modules="numpy")
        
        def f(q):
            return lam(*q)
        return f