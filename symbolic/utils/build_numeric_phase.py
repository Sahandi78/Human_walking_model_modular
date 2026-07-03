import sympy as sp
import numpy as np

# Utilities to lambdify symbolic expressions into fast numeric functions

from symbolic.utils.lambdify_helpers import (
    lambdify_position,
    lambdify_velocity, 
    lambdify_matrix
)
# ------------------------------------------------------------
# Build numeric (callable) functions for each phase
#
#    This converts symbolic expressions like:
#       XA(q), YA(q), VA(q, dq)
#
#    into fast NumPy-callable functions:
#       XA(t, y), YA(t, y), VA(t, y)
#
#    IMPORTANT:
#    - Works for both single- and double-stance phases
#    - Uses phase["q"] and phase["dq"] to infer state layout
# ------------------------------------------------------------

def build_numeric_phase(phase_sub):
    numeric = {
        "name": phase_sub["name"],
        "positions": {},
        "velocities": {},
        "dynamics": {},
        "contact": {},
        "constraint": {}
    }

    # ---- Positions ----
    for key, (X, Y) in phase_sub.get("positions", {}).items():
        numeric["positions"][key] = (
            lambdify_position(X, phase_sub),
            lambdify_position(Y, phase_sub),
        )

    # ---- Velocities ----
    for key, V in phase_sub.get("velocities", {}).items():
        numeric["velocities"][key] = lambdify_velocity(V, phase_sub)

    # ---- Dynamics ----
    numeric["dynamics"]["M"] = lambdify_matrix(
        phase_sub["M"], phase_sub, use_dq=False
    )

    numeric["dynamics"]["C"] = lambdify_matrix(
        phase_sub["C"], phase_sub, use_dq=True
    )

    numeric["dynamics"]["G"] = lambdify_matrix(
        phase_sub["G"], phase_sub, use_dq=False
    )
    
    numeric["dynamics"]["F_generalized"] = lambdify_matrix(
        phase_sub["F_generalized"], phase_sub, use_dq=True
    )

    # ---- Constraint Jacobian (double-stance only) ----
    if phase_sub.get("D") is not None:
        numeric["dynamics"]["D"] = lambdify_matrix(
            phase_sub["D"], phase_sub, use_dq=False
        )

    # ---- Constraint acceleration term (double-stance only) ----
    if phase_sub.get("Gamma") is not None:
        numeric["dynamics"]["Gamma"] = lambdify_matrix(
            phase_sub["Gamma"], phase_sub, use_dq=True
        )

    # ---- Contact forces ----
    for key, F in phase_sub.get("contact", {}).items():
        numeric["contact"][key] = lambdify_matrix(
            F, phase_sub, use_dq=True
        )
        
    # ---- Constraint 1 function (double-stance only) ----
    if phase_sub.get("Phi_1") is not None:
        numeric["constraint"]["Phi_1"] = lambdify_matrix(
            phase_sub["Phi_1"], phase_sub, use_dq=False
        )
        
    # ---- Constraint 2 function (double-stance only) ----
    if phase_sub.get("Phi_2") is not None:
        numeric["constraint"]["Phi_2"] = lambdify_matrix(
            phase_sub["Phi_2"], phase_sub, use_dq=False
        )

    return numeric
