"""
Gait phase transition events (numeric only)
"""

from .thresholds import (
    GROUND_HEIGHT,
    H_TOL,
    V_STRIKE,
    V_OFF,
    FORCE_OFF
)


# ------------------------------------------------------------
# STRIKE EVENTS
# ------------------------------------------------------------

def left_heel_strike(Y_lheel, dY_lheel):
    return (
        Y_lheel < GROUND_HEIGHT + H_TOL and
        dY_lheel < V_STRIKE
    )


def right_heel_strike(Y_rheel, dY_rheel):
    return (
        Y_rheel < GROUND_HEIGHT + H_TOL and
        dY_rheel < V_STRIKE
    )


def left_toe_strike(Y_lt, dY_lt):
    return (
        Y_lt < GROUND_HEIGHT + H_TOL and
        dY_lt < V_STRIKE
    )


def right_toe_strike(Y_rt, dY_rt):
    return (
        Y_rt < GROUND_HEIGHT + H_TOL and
        dY_rt < V_STRIKE
    )


# ------------------------------------------------------------
# OFF EVENTS
# ------------------------------------------------------------

def left_toe_off(Y_lt, dY_lt, F_lt=None):
    height_ok = Y_lt > GROUND_HEIGHT - H_TOL
    velocity_ok = dY_lt > V_OFF
    force_ok = (F_lt is None) or (F_lt < FORCE_OFF)

    return height_ok and velocity_ok and force_ok


def right_toe_off(Y_rt, dY_rt, F_rt=None):
    height_ok = Y_rt > GROUND_HEIGHT - H_TOL
    velocity_ok = dY_rt > V_OFF
    force_ok = (F_rt is None) or (F_rt < FORCE_OFF)

    return height_ok and velocity_ok and force_ok


def left_heel_off(Y_lheel, dY_lheel, F_lh=None):
    height_ok = Y_lheel > GROUND_HEIGHT - H_TOL
    velocity_ok = dY_lheel > V_OFF
    force_ok = (F_lh is None) or (F_lh < FORCE_OFF)

    return height_ok and velocity_ok and force_ok


def right_heel_off(Y_rheel, dY_rheel, F_rh=None):
    height_ok = Y_rheel > GROUND_HEIGHT - H_TOL
    velocity_ok = dY_rheel > V_OFF
    force_ok = (F_rh is None) or (F_rh < FORCE_OFF)

    return height_ok and velocity_ok and force_ok
