import sympy as sp

def substitute_phase(phase_dict, values):
    """
    Substitute numerical values into a phase dictionary.
    Works for single- and double-stance phases.
    """

    def safe_sub(expr):
        if expr is None:
            return None
        if isinstance(expr, (sp.Matrix, sp.Expr)):
            return expr.subs(values)
        return expr

    phase_sub = {
        "name": phase_dict["name"],

        # Generalized coordinates (unchanged symbols)
        "q": phase_dict["q"],
        "dq": phase_dict["dq"],
        "ddq": phase_dict["ddq"],

        # Dynamics
        "M": safe_sub(phase_dict.get("M")),
        "C": safe_sub(phase_dict.get("C")),
        "G": safe_sub(phase_dict.get("G")),
        "D": safe_sub(phase_dict.get("D")),
        "Gamma": safe_sub(phase_dict.get("Gamma")),
        "F_generalized": safe_sub(phase_dict.get("F_generalized")),
    }
    
    if phase_dict.get("Phi_1") is not None:
        phase_sub["Phi_1"] = safe_sub(phase_dict.get("Phi_1"))
        phase_sub["Phi_2"] = safe_sub(phase_dict.get("Phi_2"))

    # Contact forces
    if "contact" in phase_dict:
        phase_sub["contact"] = {
            k: safe_sub(v) for k, v in phase_dict["contact"].items()
        }

    # Positions
    if "positions" in phase_dict:
        phase_sub["positions"] = {
            k: (safe_sub(v[0]), safe_sub(v[1]))
            for k, v in phase_dict["positions"].items()
        }

    # Velocities
    if "velocities" in phase_dict:
        phase_sub["velocities"] = {
            k: safe_sub(v) for k, v in phase_dict["velocities"].items()
        }

    return phase_sub
