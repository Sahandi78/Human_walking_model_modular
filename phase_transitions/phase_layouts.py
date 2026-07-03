# phase_layouts.py

PHASE_LAYOUTS = {
    "S_L_T": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["L_toe"]
    },
    "S_L_H": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["L_heel"]
    },
    "S_L_TH": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["L_toe"]  # heel+toe constraint already internal
    },
    "S_R_T": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["R_toe"]
    },
    "S_R_H": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["R_heel"]
    },
    "S_R_TH": {
        "n_feet": 1,
        "has_lambda": False,
        "foot_labels": ["R_toe"]
    },
    "D_R_H": {
        "n_feet": 2,
        "has_lambda": True,
        "foot_labels": ["L_toe", "R_heel"]
    },
    "D_R_TH": {
        "n_feet": 2,
        "has_lambda": True,
        "foot_labels": ["L_toe", "R_toe"]
    },
    "D_L_H": {
        "n_feet": 2,
        "has_lambda": True,
        "foot_labels": ["R_toe", "L_heel"]
    },
    "D_L_TH": {
        "n_feet": 2,
        "has_lambda": True,
        "foot_labels": ["R_toe", "L_toe"]
    }
}
