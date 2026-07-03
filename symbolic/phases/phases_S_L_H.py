import sympy as sp
from sympy import Matrix, zeros, sin, cos, diff, symbols


#-------- Left stance human (toe off ground, heel on ground) --------
def derive():
    """
    Derives symbolic equations for phase S_L_T
    Returns:
        dict with keys: M, C, G, F, q, dq, ddq
    """

    # Defining symbolic variables
    q1, q2, q3, q4, q5, q6, q7 = sp.symbols('q1 q2 q3 q4 q5 q6 q7')
    dq1, dq2, dq3, dq4, dq5, dq6, dq7 = sp.symbols('dq1 dq2 dq3 dq4 dq5 dq6 dq7')
    ddq1, ddq2, ddq3, ddq4, ddq5, ddq6, ddq7 = sp.symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7')
    X_lheel, Y_lheel, dX_lheel, dY_lheel, ddX_lheel, ddY_lheel = sp.symbols('X_lheel Y_lheel dX_lheel dY_lheel ddX_lheel ddY_lheel')

    # Defining torques
    tau_lheel, tau_la, tau_lk, tau_lh, tau_rt, tau_ra, tau_rk, tau_rh = sp.symbols('tau_lheel tau_la tau_lk tau_lh tau_rt tau_ra tau_rk tau_rh')

    # Defining physical parameters
    m_f, m_s, m_t, m_h = sp.symbols('m_f m_s m_t m_h')
    L_t, L_f, L_s, L_h = sp.symbols('L_t L_f L_s L_h')
    r_s, r_f, r_t, r_h = sp.symbols('r_s r_f r_t r_h')
    g = sp.symbols('g')

    # Defining inertia and other parameters
    I_f, I_s, I_t, I_h = sp.symbols('I_f I_s I_t I_h')
    alpha, beta, gama, L_heelTotoe = sp.symbols('alpha beta gama L_heelTotoe')

    # Defining spring and damper constants
    k_lt, k_la, k_lk, k_lh = sp.symbols('k_lt k_la k_lk k_lh')
    K_rt, k_ra, k_rk, k_rh = sp.symbols('k_rt k_ra k_rk k_rh')
    C_lt, C_la, C_lk, C_lh = sp.symbols('C_lt C_la C_lk C_lh')
    C_rt, C_ra, C_rk, C_rh = sp.symbols('C_rt C_ra C_rk C_rh')

    # Parameters for Jacobian
    rs_f, rs_s, rs_t, rs_h = sp.symbols('rs_f rs_s rs_t rs_h')

    # Additional parameters
    h = sp.symbols('h')

    # Hunt-Crossley based normal force parameters
    k_hc = sp.Symbol('k_hc')      # Hunt-Crossley stiffness
    c_hc = sp.Symbol('c_hc')      # Hunt-Crossley damping coefficient
    n_hc = sp.Symbol('n_hc')      # Nonlinear exponent (typically 1.5 for Hertzian contact)
    eps_v = sp.symbols('eps_v')
    eps_phi = sp.symbols('eps_phi')
    v_eps = sp.symbols('v_eps')


    # friction model (LuGre-inspired) parameters
    mu_s = sp.Symbol('mu_s')  # Static friction coefficient
    mu_d = sp.Symbol('mu_d')  # Dynamic friction coefficient
    v_s = sp.Symbol('v_s')    # Stribeck velocity

    #------- Left foot calculations -------
    V_lheel__S_L_H = sp.Matrix([dX_lheel, dY_lheel, 0])
    a_lheel__S_L_H = sp.Matrix([ddX_lheel, ddY_lheel, 0])
    Omega_lf = sp.Matrix([0, 0, dq1])
    alpha_lf = sp.Matrix([0, 0, ddq1])

    # Relative positions
    X_lheel__S_L_H = X_lheel
    Y_lheel__S_L_H = Y_lheel
    XA__S_L_H = X_lheel__S_L_H + L_heelTotoe * sp.cos(sp.pi/2 - alpha + q1)
    YA__S_L_H = Y_lheel__S_L_H + L_heelTotoe * sp.sin(sp.pi/2 - alpha + q1)
    rA__S_L_H = sp.Matrix([XA__S_L_H, YA__S_L_H, 0])
    X_lbone__S_L_H = XA__S_L_H
    Y_lbone__S_L_H = YA__S_L_H
    V_lbone__S_L_H = V_lheel__S_L_H + Omega_lf.cross(sp.Matrix([L_heelTotoe * sp.cos(sp.pi/2 - alpha + q1), L_heelTotoe * sp.sin(sp.pi/2 - alpha + q1), 0]))
    VA__S_L_H = V_lbone__S_L_H

    # RElative positions
    rG_lf = sp.Matrix([XA__S_L_H - r_f * sp.sin(q1), YA__S_L_H + r_f * sp.cos(q1), 0])
    rBA = sp.Matrix([-L_f * sp.sin(q1), L_f * sp.cos(q1), 0])
    rIA = sp.Transpose(sp.Matrix([-rs_f * sp.sin(q1), rs_f * sp.cos(q1), 0]))

    # Velocity and acceleration
    VG_lf = V_lheel__S_L_H + Omega_lf.cross(rG_lf)
    VB = V_lheel__S_L_H + Omega_lf.cross(rBA)
    aG_lf = a_lheel__S_L_H + alpha_lf.cross(rG_lf) + Omega_lf.cross(Omega_lf.cross(rG_lf))
    aA = a_lheel__S_L_H + alpha_lf.cross(rA__S_L_H) + Omega_lf.cross(Omega_lf.cross(rA__S_L_H))
    aB = aA + alpha_lf.cross(rBA) + Omega_lf.cross(Omega_lf.cross(rBA))

    # Energy calculations
    T1 = 1/2 * m_f * (VG_lf.dot(VG_lf)) + 1/2 * I_f * (Omega_lf.dot(Omega_lf))
    V1 = m_f * g * r_f * sp.cos(q1)

    # GRF for left heel
    # Penetration depth and rate (compression positive)
    delta_heel     = sp.Max(-Y_lheel, 0)
    contact_active = (1 + sp.tanh(delta_heel / eps_phi)) / 2
    vel_smooth = sp.tanh((-dY_lheel) / v_eps)
    delta_dot_heel = (-dY_lheel) * contact_active * vel_smooth
    # Hunt-Crossley normal force
    # F_n = k * delta^n * (1 + c * delta_dot)
    # The damping term is proportional to penetration depth (key feature!)
    Fy_heel_hc = k_hc * delta_heel**n_hc * (1 + c_hc * delta_dot_heel)
    # Stribeck friction curve for smoother transition
    mu_eff = mu_d + (mu_s - mu_d) * sp.exp(-(dX_lheel/v_s)**2)
    # Friction force with smooth saturation
    Fx_heel_hc = -mu_eff * Fy_heel_hc * sp.tanh(dX_lheel / eps_v)
    # Combined ground reaction force vector
    F_heel__S_L_H = sp.Matrix([Fx_heel_hc, Fy_heel_hc, 0])

    #-------- Left shank calculations --------
    Omega_ls = sp.Matrix([0, 0, dq2])
    alpha_ls = sp.Matrix([0, 0, ddq2])

    # Relative positions for left shank
    rG_ls = sp.Matrix([-r_s * sp.sin(q2), r_s * sp.cos(q2), 0])
    rCB = sp.Matrix([-L_s * sp.sin(q2), L_s * sp.cos(q2), 0])
    rJB = sp.Transpose(sp.Matrix([-rs_s * sp.sin(q2), rs_s * sp.cos(q2), 0]))

    # Velocity and acceleration for left shank
    VG_ls = VB + Omega_ls.cross(rG_ls)
    VC = VB + Omega_ls.cross(rCB)
    aG_ls = aB + alpha_ls.cross(rG_ls) + Omega_ls.cross(Omega_ls.cross(rG_ls))
    aC = aB + alpha_ls.cross(rCB) + Omega_ls.cross(Omega_ls.cross(rCB))

    # Energy calculations for left shank
    T2 = 1/2 * m_s * (VG_ls.dot(VG_ls)) + 1/2 * I_s * (Omega_ls.dot(Omega_ls))
    V2 = m_s * g * (Y_lheel + L_f * sp.cos(q1) + r_s * sp.cos(q2))
    Ve_la = 1/2 * k_la * (q2 - q1)**2
    R_la = 1/2 * C_la * ((Omega_ls - Omega_lf).dot(Omega_ls - Omega_lf))

    #-------- Left thigh calculations --------
    Omega_lt = sp.Matrix([0, 0, dq3])
    alpha_lt = sp.Matrix([0, 0, ddq3])

    # Relative positions for left thigh
    rG_lt = sp.Matrix([-r_t * sp.sin(q3), r_t * sp.cos(q3), 0])
    rDC = sp.Matrix([-L_t * sp.sin(q3), L_t * sp.cos(q3), 0])
    rKC = sp.Transpose(sp.Matrix([-rs_t * sp.sin(q3), rs_t * sp.cos(q3), 0]))

    # Velocity and acceleration for left thigh
    VG_lt = VC + Omega_lt.cross(rG_lt)
    VD = VC + Omega_lt.cross(rDC)
    aG_lt = aC + alpha_lt.cross(rG_lt) + Omega_lt.cross(Omega_lt.cross(rG_lt))
    aD = aC + alpha_lt.cross(rDC) + Omega_lt.cross(Omega_lt.cross(rDC))

    # Energy calculations for left thigh
    T3 = 1/2 * m_t * (VG_lt.dot(VG_lt)) + 1/2 * I_t * (Omega_lt.dot(Omega_lt))
    V3 = m_t * g * (Y_lheel + L_f * sp.cos(q1) + L_s * sp.cos(q2) + r_t * sp.cos(q3))
    Ve_lk = 1/2 * k_lk * (q3 - q2)**2
    R_lk = 1/2 * C_lk * ((Omega_lt - Omega_ls).dot(Omega_lt - Omega_ls))

    #-------- HAT (Head, Arms, Torso) --------
    Omega_h = sp.Matrix([0, 0, dq4])
    alpha_h = sp.Matrix([0, 0, ddq4])

    # Relative positions for HAT
    rG_h = sp.Matrix([-r_h * sp.sin(q4), r_h * sp.cos(q4), 0])
    rED = sp.Matrix([-L_h * sp.sin(q4), L_h * sp.cos(q4), 0])
    rLD = sp.Transpose(sp.Matrix([-rs_h * sp.sin(q4), rs_h * sp.cos(q4), 0]))

    # Velocity and acceleration for HAT
    VG_h = VD + Omega_h.cross(rG_h)
    VE = VD + Omega_h.cross(rED)
    aG_h = aD + alpha_h.cross(rG_h) + Omega_h.cross(Omega_h.cross(rG_h))
    aE = aD + alpha_h.cross(rED) + Omega_h.cross(Omega_h.cross(rED))

    # Energy calculations for HAT
    T4 = 1/2 * m_h * (VG_h.dot(VG_h)) + 1/2 * I_h * (Omega_h.dot(Omega_h))
    V4 = m_h * g * (Y_lheel + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) + r_h * sp.cos(q4))
    Ve_lh = 1/2 * k_lh * (q4 - q3)**2
    R_lh = 1/2 * C_lh * ((Omega_h - Omega_lt).dot(Omega_h - Omega_lt))

    #-------- Right thigh calculations --------
    s_t = L_t - r_t
    ss_t = L_t - rs_t
    Omega_rt = sp.Matrix([0, 0, dq5])
    alpha_rt = sp.Matrix([0, 0, ddq5])

    # Relative positions for right thigh
    rG_rt = sp.Matrix([s_t * sp.sin(q5), -s_t * sp.cos(q5), 0])
    rFD = sp.Matrix([L_t * sp.sin(q5), -L_t * sp.cos(q5), 0])
    rMD = sp.Transpose(sp.Matrix([ss_t * sp.sin(q5), -ss_t * sp.cos(q5), 0]))

    # Velocity and acceleration for right thigh
    VG_rt = VD + Omega_rt.cross(rG_rt)
    VF = VD + Omega_rt.cross(rFD)
    aG_rt = aD + alpha_rt.cross(rG_rt) + Omega_rt.cross(Omega_rt.cross(rG_rt))
    aF = aD + alpha_rt.cross(rFD) + Omega_rt.cross(Omega_rt.cross(rFD))

    # Energy calculations for right thigh
    T5 = 1/2 * m_t * (VG_rt.dot(VG_rt)) + 1/2 * I_t * (Omega_rt.dot(Omega_rt))
    V5 = m_t * g * (Y_lheel + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) - s_t * sp.cos(q5))
    Ve_rh = 1/2 * k_rh * (q5 - q4)**2
    R_rh = 1/2 * C_rh * ((Omega_rt - Omega_h).dot(Omega_rt - Omega_h))

    #-------- Right shank calculations --------
    s_s = L_s - r_s
    ss_s = L_s - rs_s
    Omega_rs = sp.Matrix([0, 0, dq6])
    alpha_rs = sp.Matrix([0, 0, ddq6])

    # Relative positions for right shank
    rG_rs = sp.Matrix([s_s * sp.sin(q6), -s_s * sp.cos(q6), 0])
    rGF = sp.Matrix([L_s * sp.sin(q6), -L_s * sp.cos(q6), 0])
    rNF = sp.Transpose(sp.Matrix([ss_s * sp.sin(q6), -ss_s * sp.cos(q6), 0]))

    # Velocity and acceleration for right shank
    VG_rs = VF + Omega_rs.cross(rG_rs)
    VG = VF + Omega_rs.cross(rGF)
    aG_rs = aF + alpha_rs.cross(rG_rs) + Omega_rs.cross(Omega_rs.cross(rG_rs))
    aG = aF + alpha_rs.cross(rGF) + Omega_rs.cross(Omega_rs.cross(rGF))

    # Energy calculations for right shank
    T6 = 1/2 * m_s * (VG_rs.dot(VG_rs)) + 1/2 * I_s * (Omega_rs.dot(Omega_rs))
    V6 = m_s * g * (Y_lheel + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) - L_t * sp.cos(q5) - s_s * sp.cos(q6))
    Ve_rk = 1/2 * k_rk * (q6 - q5)**2
    R_rk = 1/2 * C_rk * ((Omega_rs - Omega_rt).dot(Omega_rs - Omega_rt))

    #-------- Right foot calculations --------
    s_f = L_f - r_f
    ss_f = L_f - rs_f
    Omega_rf = sp.Matrix([0, 0, dq7])
    alpha_rf = sp.Matrix([0, 0, ddq7])

    # Relative positions for right foot
    rG_rf = sp.Matrix([s_f * sp.sin(q7), -s_f * sp.cos(q7), 0])
    rHG = sp.Matrix([L_f * sp.sin(q7), -L_f * sp.cos(q7), 0])
    rOG = sp.Transpose(sp.Matrix([ss_f * sp.sin(q7), -ss_f * sp.cos(q7), 0]))

    # Velocity and acceleration for right foot
    VG_rf = VG + Omega_rf.cross(rG_rf)
    VH__S_L_H = VG + Omega_rf.cross(rHG)
    aG_rf = aG + alpha_rf.cross(rG_rf) + Omega_rf.cross(Omega_rf.cross(rG_rf))
    aH = aG + alpha_rf.cross(rHG) + Omega_rf.cross(Omega_rf.cross(rHG))

    # Energy calculations for right foot
    T7 = 1/2 * m_f * (VG_rf.dot(VG_rf)) + 1/2 * I_f * (Omega_rf.dot(Omega_rf))
    V7 = m_f * g * (Y_lheel + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) - L_t * sp.cos(q5) - L_s * sp.cos(q6) - s_f * sp.cos(q7))
    Ve_ra = 1/2 * k_ra * (q7 - q6)**2
    R_ra = 1/2 * C_ra * ((Omega_rf - Omega_rs).dot(Omega_rf - Omega_rs))


    # Positions of points B, C, D, E, F, G, H, and right heel
    XB__S_L_H = XA__S_L_H - L_f * sp.sin(q1)
    YB__S_L_H = YA__S_L_H + L_f * sp.cos(q1)
    rB__S_L_H = sp.Transpose(sp.Matrix([XB__S_L_H, YB__S_L_H, 0]))

    XC__S_L_H = XB__S_L_H - L_s * sp.sin(q2)
    YC__S_L_H = YB__S_L_H + L_s * sp.cos(q2)
    rC__S_L_H = sp.Transpose(sp.Matrix([XC__S_L_H, YC__S_L_H, 0]))

    XD__S_L_H = XC__S_L_H - L_t * sp.sin(q3)
    YD__S_L_H = YC__S_L_H + L_t * sp.cos(q3)
    rD__S_L_H = sp.Transpose(sp.Matrix([XD__S_L_H, YD__S_L_H, 0]))

    XE__S_L_H = XD__S_L_H - L_h * sp.sin(q4)
    YE__S_L_H = YD__S_L_H + L_h * sp.cos(q4)
    rE__S_L_H = sp.Transpose(sp.Matrix([XE__S_L_H, YE__S_L_H, 0]))

    XF__S_L_H = XD__S_L_H + L_t * sp.sin(q5)
    YF__S_L_H = YD__S_L_H - L_t * sp.cos(q5)
    rF__S_L_H = sp.Transpose(sp.Matrix([XF__S_L_H, YF__S_L_H, 0]))

    XG__S_L_H = XF__S_L_H + L_s * sp.sin(q6)
    YG__S_L_H = YF__S_L_H - L_s * sp.cos(q6)
    rG__S_L_H = sp.Transpose(sp.Matrix([XG__S_L_H, YG__S_L_H, 0]))

    XH__S_L_H = XG__S_L_H + L_f * sp.sin(q7)
    YH__S_L_H = YG__S_L_H - L_f * sp.cos(q7)
    rH__S_L_H = sp.Transpose(sp.Matrix([XH__S_L_H, YH__S_L_H, 0]))

    X_rheel__S_L_H = XH__S_L_H + L_heelTotoe * sp.sin(q7 - alpha)
    Y_rheel__S_L_H = YH__S_L_H - L_heelTotoe * sp.cos(q7 - alpha)

    X_rbone__S_L_H = XH__S_L_H
    Y_rbone__S_L_H = YH__S_L_H

    V_rheel__S_L_H = VH__S_L_H + Omega_rf.cross(sp.Matrix([L_heelTotoe * sp.sin(q7 - alpha), - L_heelTotoe * sp.cos(q7 - alpha), 0]))


    # Define symbols for generalized coordinates and their derivatives
    q = list(symbols('q1 q2 q3 q4 q5 q6 q7 X_lheel Y_lheel'))
    dq = list(symbols('dq1 dq2 dq3 dq4 dq5 dq6 dq7 dX_lheel dY_lheel'))
    ddq = list(symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 ddX_lheel ddY_lheel'))

    # Lagrange’s Equations of Motion
    T = T1 + T2 + T3 + T4 + T5 + T6 + T7
    Ve = Ve_la + Ve_lk + Ve_lh + Ve_rh + Ve_rk + Ve_ra
    V = V1 + V2 + V3 + V4 + V5 + V6 + V7 + Ve

    L = T - V

    # Derivatives of Lagrangian with respect to generalized velocities
    dL_dq = [sp.diff(L, dq_i) for dq_i in dq]

    # Mass matrix M
    M = sp.Matrix([[sp.diff(dL_dq_i, dq_j) for dq_j in dq] for dL_dq_i in dL_dq])

    # Coriolis and centrifugal matrix N
    C = zeros(len(q), len(q))  # Initialize Coriolis matrix
    for i in range(len(q)):
        for j in range(len(q)):
            for k in range(len(q)):
                C[i, j] += (diff(M[i, j], q[k]) +
                            diff(M[i, k], q[j]) -
                            diff(M[j, k], q[i])) * dq[k] / 2

    # Compute gradient of potential energy
    G = sp.Matrix([sp.diff(V, q_i) for q_i in q])

    # External torques
    # Work done by external forces
    dW = tau_lheel * Omega_lf[2] + tau_la * (Omega_ls[2] - Omega_lf[2]) + tau_lk * (Omega_lt[2] - Omega_ls[2]) + \
         tau_lh * (Omega_h[2] - Omega_lt[2]) + tau_rh * (Omega_rt[2] - Omega_h[2]) + tau_rk * (Omega_rs[2] - Omega_rt[2]) + \
         tau_ra * (Omega_rf[2] - Omega_rs[2])

    # Rayleigh dissipation function
    R = R_la + R_lk + R_lh + R_rh + R_rk + R_ra

    # Generalized forces Q
    Q = sp.Matrix([sp.diff(dW, dq_i) for dq_i in dq])

    # Jacobian of heel position with respect to generalized coordinates
    J_heel = sp.Matrix([X_lheel__S_L_H, Y_lheel__S_L_H, 0]).jacobian(sp.Matrix(q))

    # Generalized force contribution from ground contact
    Q_ground = J_heel.T * F_heel__S_L_H

    # Total generalized forces
    Q_total_ls = Q + Q_ground

    # Velocity-dependent generalized forces Q_v
    Q_v = -sp.Matrix([sp.diff(R, dq_i) for dq_i in dq])

    # Damping matrix C_v
    C_v = sp.Matrix([[sp.diff(Q_v_i, dq_j) for dq_j in dq] for Q_v_i in Q_v])

    # Total energy E
    E_ls = T + V

    # Time derivative of energy dE_dt
    dE_dt_ls = sum(sp.Matrix([sp.diff(E_ls, q_i) * dq_i + sp.diff(E_ls, dq_i) * ddq_i for q_i, dq_i, ddq_i in zip(q, dq, ddq)]))

    # Work done by non-conservative forces dW_nc
    dW_nc = sum([sp.diff(Q_i, dq_j) * dq_j for Q_i, dq_j in zip(Q, dq)]) - sum([sp.diff(R, dq_i) * dq_i for dq_i in dq])


    return {
    "name": "S_L_H",

    # Generalized coordinates
    "q": q,
    "dq": dq,
    "ddq": ddq,

    # Dynamics
    "M": M,
    "C": C + C_v,
    "G": G,
    "F_generalized": Q_ground,

    # Contact forces (raw)
    "contact": {
        "heel_force": F_heel__S_L_H
    },

    # Positions (symbolic, BEFORE substitution)
    "positions": {
        "A": (XA__S_L_H, YA__S_L_H),
        "B": (XB__S_L_H, YB__S_L_H),
        "C": (XC__S_L_H, YC__S_L_H),
        "D": (XD__S_L_H, YD__S_L_H),
        "E": (XE__S_L_H, YE__S_L_H),
        "F": (XF__S_L_H, YF__S_L_H),
        "G": (XG__S_L_H, YG__S_L_H),
        "H": (XH__S_L_H, YH__S_L_H),
        "lheel": (X_lheel__S_L_H, Y_lheel__S_L_H),
        "rheel": (X_rheel__S_L_H, Y_rheel__S_L_H)
    },

    # Velocities
    "velocities": {
        "A": VA__S_L_H,
        "lheel": V_lheel__S_L_H,
        "H": VH__S_L_H,
        "rheel": V_rheel__S_L_H
    }
}
