from events.gait_events import *
from phase_transitions.transitions import transition_phase
from simulation.grf import compute_grf
import numpy as np

def update_phase(t, y_current, current_phase, numeric_phases, F_lt, F_lh, F_rt, F_rh):
    """
    Detect events and apply phase transitions
    Returns: y_current (possibly updated), current_phase
    """
    F_lt_y = float(F_lt[1])
    F_lh_y = float(F_lh[1])
    F_rt_y = float(F_rt[1])
    F_rh_y = float(F_rh[1])
    
    numeric_phase = numeric_phases[current_phase]
    
    XA_fun, YA_fun = numeric_phase["positions"]["A"]
    XA = XA_fun(t, y_current)
    YA = YA_fun(t, y_current)
    
    XH_fun, YH_fun = numeric_phase["positions"]["H"]
    XH = XH_fun(t, y_current)
    YH = YH_fun(t, y_current)

    X_lheel_fun, Y_lheel_fun = numeric_phase["positions"]["lheel"]
    X_lheel = X_lheel_fun(t, y_current)
    Y_lheel = Y_lheel_fun(t, y_current)
    
    X_rheel_fun, Y_rheel_fun = numeric_phase["positions"]["rheel"]
    X_rheel = X_rheel_fun(t, y_current)
    Y_rheel = Y_rheel_fun(t, y_current)

    
    VA = numeric_phase["velocities"]["A"](t, y_current)
    VH = numeric_phase["velocities"]["H"](t, y_current)
    V_lheel = numeric_phase["velocities"]["lheel"](t, y_current)
    V_rheel = numeric_phase["velocities"]["rheel"](t, y_current)

    # S_L_TH → S_L_T
    if current_phase == 'S_L_TH' and left_heel_off(Y_lheel, V_lheel[1], F_lh_y):
        print(f"Left heel off detected at t = {t}")
        current_phase = 'S_L_T'
        print(current_phase)
        print('========================================')
    
    # S_L_T → D_R_H
    elif current_phase == 'S_L_T' and right_heel_strike(Y_rheel, V_rheel[1]):
        print(f"Right heel strike detected at t = {t}")
        current_phase = 'D_R_H'
        y_current = transition_phase(
        t,
        y_current,
        'S_L_T',
        'D_R_H',
        numeric_phases,
        new_contacts={
        'L_toe': get_point_state('A', t, y_current, numeric_phase),
        'R_heel': get_point_state('rheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # D_R_H → D_R_TH
    elif current_phase == 'D_R_H' and right_toe_strike(YH, VH[1]):
        print(f"Right toe strike detected at t = {t}")
        last_phase = current_phase
        current_phase = 'D_R_TH'
        print(current_phase)
        print('========================================')
   
    # D_R_H → S_R_H
    elif current_phase == 'D_R_H' and left_toe_off(YA, VA[1], F_lt_y):
        print(f"Left toe off detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_R_H'
        y_current = transition_phase(
        t, y_current,
        'D_R_H',
        'S_R_H',
        numeric_phases,
        new_contacts={
        'R_heel': get_point_state('rheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # S_R_H → S_R_TH
    elif current_phase == 'S_R_H' and right_toe_strike(YH, VH[1]):
        print(f"Right toe strike detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_R_TH'
        y_current = transition_phase(
        t,
        y_current,
        "S_R_H",
        "S_R_TH",
        numeric_phases,
        new_contacts={
        'R_toe': get_point_state('H', t, y_current, numeric_phase),
        'R_heel': get_point_state('rheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # D_R_TH → S_R_TH
    elif current_phase == 'D_R_TH' and left_toe_off(YA, VA[1], F_lt_y):
        print(f"Left toe off detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_R_TH'  
        y_current = transition_phase(
        t,
        y_current,
        "D_R_TH",
        "S_R_TH",
        numeric_phases,
        new_contacts={
        'R_toe': get_point_state('H', t, y_current, numeric_phase),
        'R_heel': get_point_state('rheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')

    # S_R_TH → S_R_T
    elif current_phase == 'S_R_TH' and right_heel_off(Y_rheel, V_rheel[1], F_rh_y):
        print(f"Right heel off detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_R_T'
        print(current_phase)
        print('========================================')
    
    # S_R_T → D_L_H
    elif current_phase == 'S_R_T' and left_heel_strike(Y_lheel, V_lheel[1]):
        print(f"Left heel strike detected at t = {t}")
        last_phase = current_phase
        current_phase = 'D_L_H'
        y_current = transition_phase(
        t,
        y_current,
        "S_R_T",
        "D_L_H",
        numeric_phases,
        new_contacts={
        'R_toe': get_point_state('H', t, y_current, numeric_phase),
        'L_heel': get_point_state('lheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # D_L_H → D_L_TH
    elif current_phase == 'D_L_H' and left_toe_strike(YA, VA[1]):
        print(f"Left toe strike detected at t = {t}")
        last_phase = current_phase
        current_phase = 'D_L_TH'
        print(current_phase)
        print('========================================')
    
    # D_L_H → S_L_H
    elif current_phase == 'D_L_H' and right_toe_off(YH, VH[1], F_rt_y):
        print(f"Right toe off detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_L_H'
        y_current = transition_phase(
        t,
        y_current,
        "D_L_H",
        "S_L_H",
        numeric_phases,
        new_contacts={
        'L_heel': get_point_state('lheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # S_L_H → S_L_TH
    elif current_phase == 'S_L_H' and left_toe_strike(YA, VA[1]):
        print(f"Left toe strike detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_L_TH'
        y_current = transition_phase(
        t,
        y_current,
        "S_L_H",
        "S_L_TH",
        numeric_phases,
        new_contacts={
        'L_toe': get_point_state('A', t, y_current, numeric_phase),
        'L_heel': get_point_state('lheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    
    # D_L_TH → S_L_TH
    elif current_phase == 'D_L_TH' and right_toe_off(YH, VH[1], F_rt_y):
        print(f"Right toe off detected at t = {t}")
        last_phase = current_phase
        current_phase = 'S_L_TH'
        y_current = transition_phase(
        t,
        y_current,
        "D_L_TH",
        "S_L_TH",
        numeric_phases,
        new_contacts={
        'L_toe': get_point_state('A', t, y_current, numeric_phase),
        'L_heel': get_point_state('rheel', t, y_current, numeric_phase)
        })
        print(current_phase)
        print('========================================')
    

    return y_current, current_phase




def get_point_state(label, t, y, numeric_phase):
    """
    Returns (X, Y, dX, dY) for a named point
    """
    X_fun, Y_fun = numeric_phase["positions"][label]
    V_fun = numeric_phase["velocities"][label]

    X = X_fun(t, y)
    Y = Y_fun(t, y)

    V = V_fun(t, y)

    # Handle scalar vs vector velocity safely
    if np.ndim(V) == 0:
        dX, dY = 0.0, V
    else:
        dX, dY = V[0].item(), V[1].item()

    return X, Y, dX, dY

