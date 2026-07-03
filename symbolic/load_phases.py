from symbolic.phases.phases_S_L_T import derive as S_L_T
from symbolic.phases.phases_S_L_H import derive as S_L_H
from symbolic.phases.phases_S_L_TH import derive as S_L_TH
from symbolic.phases.phases_D_R_H import derive as D_R_H
from symbolic.phases.phases_D_R_TH import derive as D_R_TH
from symbolic.phases.phases_S_R_T import derive as S_R_T
from symbolic.phases.phases_S_R_H import derive as S_R_H
from symbolic.phases.phases_S_R_TH import derive as S_R_TH
from symbolic.phases.phases_D_L_H import derive as D_L_H
from symbolic.phases.phases_D_L_TH import derive as D_L_TH


def load_all_phases():
    phases = {}

    for derive_fn in [
        S_L_T,
        S_L_H,
        S_L_TH,
        D_R_H, 
        D_R_TH,
        S_R_T,
        S_R_H,
        S_R_TH,
        D_L_H, 
        D_L_TH
    ]:
        phase = derive_fn()
        phases[phase["name"]] = phase

    return phases
