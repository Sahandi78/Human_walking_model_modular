import sympy as sp
from sympy import Matrix, zeros, sin, cos, diff, symbols


#-------Right stance human (toe on ground, heel off ground)-------
def derive():
    """
    Derives symbolic equations for phase S_L_T
    Returns:
        dict with keys: M, C, G, F, q, dq, ddq
    """

    # Define symbols
    q1, q2, q3, q4, q5, q6, q7 = sp.symbols('q1 q2 q3 q4 q5 q6 q7')
    dq1, dq2, dq3, dq4, dq5, dq6, dq7 = sp.symbols('dq1 dq2 dq3 dq4 dq5 dq6 dq7')
    ddq1, ddq2, ddq3, ddq4, ddq5, ddq6, ddq7 = sp.symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7')
    X_rt, Y_rt, dX_rt, dY_rt, ddX_rt, ddY_rt = sp.symbols('X_rt Y_rt dX_rt dY_rt ddX_rt ddY_rt')

    # Define torques
    tau_lt, tau_la, tau_lk, tau_lh, tau_rt, tau_ra, tau_rk, tau_rh = sp.symbols('tau_lt tau_la tau_lk tau_lh tau_rt tau_ra tau_rk tau_rh')

    # Define masses, lengths, distances to gravity center, and gravity
    m_to, m_f, m_s, m_t, m_h = sp.symbols('m_to m_f m_s m_t m_h')
    L_to, L_t, L_f, L_s, L_h = sp.symbols('L_to L_t L_f L_s L_h')
    r_to, r_s, r_f, r_t, r_h = sp.symbols('r_to r_s r_f r_t r_h')
    g = sp.symbols('g')

    # Define inertias and foot parameters
    I_f, I_s, I_t, I_h = sp.symbols('I_f I_s I_t I_h')
    alpha, beta, gama, L_heelTotoe = sp.symbols('alpha beta gama L_heelTotoe')

    # Define spring and damper constants
    k_lt, k_la, k_lk, k_lh = sp.symbols('k_lt k_la k_lk k_lh')
    k_rt, k_ra, k_rk, k_rh = sp.symbols('k_rt k_ra k_rk k_rh')
    C_lt, C_la, C_lk, C_lh = sp.symbols('C_lt C_la C_lk C_lh')
    C_rt, C_ra, C_rk, C_rh = sp.symbols('C_rt C_ra C_rk C_rh')

    # Define Jacobian parameters
    rs_f, rs_s, rs_t, rs_h = sp.symbols('rs_f rs_s rs_t rs_h')

    # Additional symbols
    h = sp.symbols('h')
    C_lheel, C_rheel, C_rbone, C_lbone = sp.symbols('C_lheel C_rheel C_rbone C_lbone')
    k_lheel, k_rheel, k_rbone, k_lbone = sp.symbols('k_lheel k_rheel k_rbone k_lbone')
    Toe_position = sp.symbols('Toe_position')

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


    # Right foot calculations
    VH__S_R_T = sp.Matrix([dX_rt, dY_rt, 0])
    aH = sp.Matrix([ddX_rt, ddY_rt, 0])
    Omega_rf = sp.Matrix([0, 0, dq7])
    alpha_rf = sp.Matrix([0, 0, ddq7])

    # Relative positions
    rG_rf = sp.Matrix([-r_f * sp.sin(q7), r_f * sp.cos(q7), 0])
    rGH = sp.Matrix([-L_f * sp.sin(q7), L_f * sp.cos(q7), 0])
    rOH = sp.Transpose(sp.Matrix([-rs_f * sp.sin(q7), rs_f * sp.cos(q7), 0]))

    VG_rf = VH__S_R_T + Omega_rf.cross(rG_rf)
    VG = VH__S_R_T + Omega_rf.cross(rGH)

    aG_rf = aH + alpha_rf.cross(rG_rf) + Omega_rf.cross(Omega_rf.cross(rG_rf))
    aG = aH + alpha_rf.cross(rGH) + Omega_rf.cross(Omega_rf.cross(rGH))

    XH__S_R_T = X_rt
    YH__S_R_T = Y_rt
    rH__S_R_T = sp.Transpose(sp.Matrix([XH__S_R_T, YH__S_R_T, 0]))
    X_rheel__S_R_T = XH__S_R_T - L_heelTotoe * sp.sin(sp.pi - alpha + q7)
    Y_rheel__S_R_T = YH__S_R_T + L_heelTotoe * sp.cos(sp.pi - alpha + q7)
    V_rheel__S_R_T = VH__S_R_T + Omega_rf.cross(sp.Matrix([-L_heelTotoe * sp.sin(sp.pi - alpha + q7), L_heelTotoe * sp.cos(sp.pi - alpha + q7), 0]))

    T7 = 1/2 * m_f * (VG_rf.dot(VG_rf)) + 1/2 * I_f * (Omega_rf.dot(Omega_rf))
    V7 = m_f * g * (Y_rt + r_f * sp.cos(q7))
    Ve_rt = 1/2 * k_rt * (alpha + q7 - sp.pi/2)**2
    R_rt = 1/2 * C_rt * (Omega_rf.dot(Omega_rf))


    # GRF for right toe
    # Penetration depth and rate (compression positive)
    delta_toe = sp.Max(-Y_rt, 0)
    contact_active = (1 + sp.tanh(delta_toe / eps_phi)) / 2
    vel_smooth = sp.tanh((-dY_rt) / v_eps)
    delta_dot_toe = (-dY_rt) * contact_active * vel_smooth
    # Hunt-Crossley normal force
    # F_n = k * delta^n * (1 + c * delta_dot)
    # The damping term is proportional to penetration depth (key feature!)
    Fy_toe_hc = k_hc * delta_toe**n_hc * (1 + c_hc * delta_dot_toe)
    # Stribeck friction curve for smoother transition
    mu_eff = mu_d + (mu_s - mu_d) * sp.exp(-(dX_rt/v_s)**2)
    # Friction force with smooth saturation
    Fx_toe_hc = -mu_eff * Fy_toe_hc * sp.tanh(dX_rt / eps_v)
    # Combined ground reaction force vector
    F_toe__S_R_T = sp.Matrix([Fx_toe_hc, Fy_toe_hc, 0])


    # Right shank calculations
    Omega_rs = sp.Matrix([0, 0, dq6])
    alpha_rs = sp.Matrix([0, 0, ddq6])

    # Relative positions
    rG_rs = sp.Matrix([-r_s * sp.sin(q6), r_s * sp.cos(q6), 0])
    rFG = sp.Matrix([-L_s * sp.sin(q6), L_s * sp.cos(q6), 0])
    rNG = sp.Transpose(sp.Matrix([-rs_s * sp.sin(q6), rs_s * sp.cos(q6), 0]))

    VG_rs = VG + Omega_rs.cross(rG_rs)
    VF = VG + Omega_rs.cross(rFG)

    aG_rs = aG + alpha_rs.cross(rG_rs) + Omega_rs.cross(Omega_rs.cross(rG_rs))
    aF = aG + alpha_rs.cross(rFG) + Omega_rs.cross(Omega_rs.cross(rFG))

    T6 = 1/2 * m_s * (VG_rs.dot(VG_rs)) + 1/2 * I_s * (Omega_rs.dot(Omega_rs))
    V6 = m_s * g * (Y_rt + L_f * sp.cos(q7) + r_s * sp.cos(q6))
    Ve_ra = 1/2 * k_ra * (q6 - q7)**2
    R_ra = 1/2 * C_ra * ((Omega_rs - Omega_rf).dot(Omega_rs - Omega_rf))

    # Right thigh calculations
    Omega_rt = sp.Matrix([0, 0, dq5])
    alpha_rt = sp.Matrix([0, 0, ddq5])

    # Relative positions
    rG_rt = sp.Matrix([-r_t * sp.sin(q5), r_t * sp.cos(q5), 0])
    rDF = sp.Matrix([-L_t * sp.sin(q5), L_t * sp.cos(q5), 0])
    rMF = sp.Transpose(sp.Matrix([-rs_t * sp.sin(q5), rs_t * sp.cos(q5), 0]))

    VG_rt = VF + Omega_rt.cross(rG_rt)
    VD = VF + Omega_rt.cross(rDF)

    aG_rt = aF + alpha_rt.cross(rG_rt) + Omega_rt.cross(Omega_rt.cross(rG_rt))
    aD = aF + alpha_rt.cross(rDF) + Omega_rt.cross(Omega_rt.cross(rDF))

    T5 = 1/2 * m_t * (VG_rt.dot(VG_rt)) + 1/2 * I_t * (Omega_rt.dot(Omega_rt))
    V5 = m_t * g * (Y_rt + L_f * sp.cos(q7) + L_s * sp.cos(q6) + r_t * sp.cos(q5))
    Ve_rk = 1/2 * k_rk * (q5 - q6)**2
    R_rk = 1/2 * C_rk * ((Omega_rt - Omega_rs).dot(Omega_rt - Omega_rs))

    # HAT calculations
    Omega_h = sp.Matrix([0, 0, dq4])
    alpha_h = sp.Matrix([0, 0, ddq4])

    # Relative positions
    rG_h = sp.Matrix([-r_h * sp.sin(q4), r_h * sp.cos(q4), 0])
    rED = sp.Matrix([-L_h * sp.sin(q4), L_h * sp.cos(q4), 0])
    rLD = sp.Transpose(sp.Matrix([-rs_h * sp.sin(q4), rs_h * sp.cos(q4), 0]))

    VG_h = VD + Omega_h.cross(rG_h)
    VE = VD + Omega_h.cross(rED)

    aG_h = aD + alpha_h.cross(rG_h) + Omega_h.cross(Omega_h.cross(rG_h))
    aE = aD + alpha_h.cross(rED) + Omega_h.cross(Omega_h.cross(rED))

    T4 = 1/2 * m_h * (VG_h.dot(VG_h)) + 1/2 * I_h * (Omega_h.dot(Omega_h))
    V4 = m_h * g * (Y_rt + L_f * sp.cos(q7) + L_s * sp.cos(q6) + L_t * sp.cos(q5) + r_h * sp.cos(q4))
    Ve_rh = 1/2 * k_rh * (q4 - q5)**2
    R_rh = 1/2 * C_rh * ((Omega_h - Omega_rt).dot(Omega_h - Omega_rt))

    # Left thigh calculations
    s_t = L_t - r_t
    ss_t = L_t - rs_t

    Omega_lt = sp.Matrix([0, 0, dq3])
    alpha_lt = sp.Matrix([0, 0, ddq3])

    # Relative positions
    rG_lt = sp.Matrix([s_t * sp.sin(q3), -s_t * sp.cos(q3), 0])
    rCD = sp.Matrix([L_t * sp.sin(q3), -L_t * sp.cos(q3), 0])
    rKD = sp.Transpose(sp.Matrix([ss_t * sp.sin(q3), -ss_t * sp.cos(q3), 0]))

    VG_lt = VD + Omega_lt.cross(rG_lt)
    VC = VD + Omega_lt.cross(rCD)

    aG_lt = aD + alpha_lt.cross(rG_lt) + Omega_lt.cross(Omega_lt.cross(rG_lt))
    aC = aD + alpha_lt.cross(rCD) + Omega_lt.cross(Omega_lt.cross(rCD))

    T3 = 1/2 * m_t * (VG_lt.dot(VG_lt)) + 1/2 * I_t * (Omega_lt.dot(Omega_lt))
    V3 = m_t * g * (Y_rt + L_f * sp.cos(q7) + L_s * sp.cos(q6) + L_t * sp.cos(q5) - s_t * sp.cos(q3))
    Ve_lh = 1/2 * k_lh * (q3 - q4)**2
    R_lh = 1/2 * C_lh * ((Omega_lt - Omega_h).dot(Omega_lt - Omega_h))

    # Left shank calculations
    s_s = L_s - r_s
    ss_s = L_s - rs_s

    Omega_ls = sp.Matrix([0, 0, dq2])
    alpha_ls = sp.Matrix([0, 0, ddq2])

    # Relative positions
    rG_ls = sp.Matrix([s_s * sp.sin(q2), -s_s * sp.cos(q2), 0])
    rBC = sp.Matrix([L_s * sp.sin(q2), -L_s * sp.cos(q2), 0])
    rJC = sp.Transpose(sp.Matrix([ss_s * sp.sin(q2), -ss_s * sp.cos(q2), 0]))

    VG_ls = VC + Omega_ls.cross(rG_ls)
    VB = VC + Omega_ls.cross(rBC)

    aG_ls = aC + alpha_ls.cross(rG_ls) + Omega_ls.cross(Omega_ls.cross(rG_ls))
    aB = aC + alpha_ls.cross(rBC) + Omega_ls.cross(Omega_ls.cross(rBC))

    T2 = 1/2 * m_s * (VG_ls.dot(VG_ls)) + 1/2 * I_s * (Omega_ls.dot(Omega_ls))
    V2 = m_s * g * (Y_rt + L_f * sp.cos(q7) + L_s * sp.cos(q6) + L_t * sp.cos(q5) - L_t * sp.cos(q3) - s_s * sp.cos(q2))
    Ve_lk = 1/2 * k_lk * (q2 - q3)**2
    R_lk = 1/2 * C_lk * ((Omega_ls - Omega_lt).dot(Omega_ls - Omega_lt))

    # Left foot calculations
    s_f = L_f - r_f
    ss_f = L_f - rs_f

    Omega_lf = sp.Matrix([0, 0, dq1])
    alpha_lf = sp.Matrix([0, 0, ddq1])

    # Relative positions
    rG_lf = sp.Matrix([s_f * sp.sin(q1), -s_f * sp.cos(q1), 0])
    rAB = sp.Matrix([L_f * sp.sin(q1), -L_f * sp.cos(q1), 0])
    rIB = sp.Transpose(sp.Matrix([ss_f * sp.sin(q1), -ss_f * sp.cos(q1), 0]))

    VG_lf = VB + Omega_lf.cross(rG_lf)
    VA__S_R_T = VB + Omega_lf.cross(rAB)

    # a = a0 + alpha * r + w * (w * r)
    aG_lf = aB + alpha_lf.cross(rG_lf) + Omega_lf.cross(Omega_lf.cross(rG_lf))
    aA = aB + alpha_lf.cross(rAB) + Omega_lf.cross(Omega_lf.cross(rAB))

    T1 = 1/2 * m_f * (VG_lf.dot(VG_lf)) + 1/2 * I_f * (Omega_lf.dot(Omega_lf))
    V1 = m_f * g * (Y_rt + L_f * sp.cos(q7) + L_s * sp.cos(q6) + L_t * sp.cos(q5) - L_t * sp.cos(q3) - L_s * sp.cos(q2) - s_f * sp.cos(q1))
    Ve_la = 1/2 * k_la * (q1 - q2)**2
    R_la = 1/2 * C_la * ((Omega_lf - Omega_ls).dot(Omega_lf - Omega_ls))

    # Position calculations
    XG__S_R_T = XH__S_R_T - L_f * sp.sin(q7)
    YG__S_R_T = YH__S_R_T + L_f * sp.cos(q7)
    rG__S_R_T = sp.Transpose(sp.Matrix([XG__S_R_T, YG__S_R_T, 0]))

    XF__S_R_T = XG__S_R_T - L_s * sp.sin(q6)
    YF__S_R_T = YG__S_R_T + L_s * sp.cos(q6)
    rF__S_R_T = sp.Transpose(sp.Matrix([XF__S_R_T, YF__S_R_T, 0]))

    XD__S_R_T = XF__S_R_T - L_t * sp.sin(q5)
    YD__S_R_T = YF__S_R_T + L_t * sp.cos(q5)
    rD__S_R_T = sp.Transpose(sp.Matrix([XD__S_R_T, YD__S_R_T, 0]))

    XE__S_R_T = XD__S_R_T - L_h * sp.sin(q4)
    YE__S_R_T = YD__S_R_T + L_h * sp.cos(q4)
    rE__S_R_T = sp.Transpose(sp.Matrix([XE__S_R_T, YE__S_R_T, 0]))

    XC__S_R_T = XD__S_R_T + L_t * sp.sin(q3)
    YC__S_R_T = YD__S_R_T - L_t * sp.cos(q3)
    rC__S_R_T = sp.Transpose(sp.Matrix([XC__S_R_T, YC__S_R_T, 0]))

    XB__S_R_T = XC__S_R_T + L_s * sp.sin(q2)
    YB__S_R_T = YC__S_R_T - L_s * sp.cos(q2)
    rB__S_R_T = sp.Transpose(sp.Matrix([XB__S_R_T, YB__S_R_T, 0]))

    XA__S_R_T = XB__S_R_T + L_f * sp.sin(q1)
    YA__S_R_T = YB__S_R_T - L_f * sp.cos(q1)
    rA__S_R_T = sp.Transpose(sp.Matrix([XA__S_R_T, YA__S_R_T, 0]))

    X_lheel__S_R_T = XA__S_R_T + L_heelTotoe * sp.sin(q1 - alpha)
    Y_lheel__S_R_T = YA__S_R_T - L_heelTotoe * sp.cos(q1 - alpha)

    X_lbone__S_R_T = XA__S_R_T
    Y_lbone__S_R_T = YA__S_R_T

    V_lheel__S_R_T = VA__S_R_T + Omega_lf.cross(sp.Matrix([L_heelTotoe * sp.sin(q1 - alpha), - L_heelTotoe * sp.cos(q1 - alpha), 0]))
    Ve_lheel = 1/2 * k_lheel * (Y_lheel__S_R_T - h)**2
    Ve_lbone = 1/2 * k_lbone * (YA__S_R_T - h)**2
    R_lheel = 1/2 * C_lheel * V_lheel__S_R_T[1]**2
    R_lbone = 1/2 * C_lbone * VA__S_R_T[1]**2

    # Define symbols for generalized coordinates and their derivatives
    q = list(symbols('q1 q2 q3 q4 q5 q6 q7 X_rt Y_rt'))
    dq = list(symbols('dq1 dq2 dq3 dq4 dq5 dq6 dq7 dX_rt dY_rt'))
    ddq = list(symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 ddX_rt ddY_rt'))

    # Lagrange’s Equations of Motion
    T = T1 + T2 + T3 + T4 + T5 + T6 + T7
    Ve = Ve_la + Ve_lk + Ve_lh + Ve_rh + Ve_rk + Ve_ra + Ve_rt
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
    dW = tau_rt * Omega_rf[2] + tau_ra * (Omega_rs[2] - Omega_rf[2]) + tau_rk * (Omega_rt[2] - Omega_rs[2]) + \
         tau_rh * (Omega_h[2] - Omega_rt[2]) + tau_lh * (Omega_lt[2] - Omega_h[2]) + tau_lk * (Omega_ls[2] - Omega_lt[2]) + \
         tau_la * (Omega_lf[2] - Omega_ls[2])

    # Rayleigh dissipation function
    R = R_la + R_lk + R_lh + R_rh + R_rk + R_ra + R_rt

    # Generalized forces Q
    Q = sp.Matrix([sp.diff(dW, dq_i) for dq_i in dq])

    # Jacobian of toe position with respect to generalized coordinates
    J_toe = sp.Matrix([XH__S_R_T, YH__S_R_T, 0]).jacobian(sp.Matrix(q))

    # Generalized force contribution from ground contact
    Q_ground = J_toe.T * F_toe__S_R_T

    # Total generalized forces
    Q_total = Q + Q_ground

    # Velocity-dependent generalized forces Q_v
    Q_v = -sp.Matrix([sp.diff(R, dq_i) for dq_i in dq])

    # Damping matrix C_v
    C_v = sp.Matrix([[sp.diff(Q_v_i, dq_j) for dq_j in dq] for Q_v_i in Q_v])

    # Total energy E
    E__S_R_T = T + V

    # Time derivative of energy dE_dt
    dE_dt_rs = sum(sp.Matrix([sp.diff(E__S_R_T, q_i) * dq_i + sp.diff(E__S_R_T, dq_i) * ddq_i for q_i, dq_i, ddq_i in zip(q, dq, ddq)]))

    # Work done by non-conservative forces dW_nc
    dWnc_dt = sum([sp.diff(Q_i, dq_j) * dq_j for Q_i, dq_j in zip(Q_total, dq)]) - sum([sp.diff(R, dq_i) * dq_i for dq_i in dq])

    # Mass Matrix difference dM
    dM = sum([diff(M, q) * dq + diff(M, dq) * ddq for q, dq, ddq in zip(q, dq, ddq)], sp.zeros(9,9))


    return {
    "name": "S_R_T",

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
        "toe_force": F_toe__S_R_T
    },

    # Positions (symbolic, BEFORE substitution)
    "positions": {
        "A": (XA__S_R_T, YA__S_R_T),
        "B": (XB__S_R_T, YB__S_R_T),
        "C": (XC__S_R_T, YC__S_R_T),
        "D": (XD__S_R_T, YD__S_R_T),
        "E": (XE__S_R_T, YE__S_R_T),
        "F": (XF__S_R_T, YF__S_R_T),
        "G": (XG__S_R_T, YG__S_R_T),
        "H": (XH__S_R_T, YH__S_R_T),
        "lheel": (X_lheel__S_R_T, Y_lheel__S_R_T),
        "rheel": (X_rheel__S_R_T, Y_rheel__S_R_T)
    },

    # Velocities
    "velocities": {
        "A": VA__S_R_T,
        "lheel": V_lheel__S_R_T,
        "H": VH__S_R_T,
        "rheel": V_rheel__S_R_T
    }
}
