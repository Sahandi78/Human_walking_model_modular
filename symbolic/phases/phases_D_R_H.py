import sympy as sp
from sympy import Matrix, zeros, sin, cos, diff, symbols


#-------Two leg stance (Left leg toe and right leg heel on ground)-------
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
    X_lt, Y_lt, dX_lt, dY_lt, ddX_lt, ddY_lt = sp.symbols('X_lt Y_lt dX_lt dY_lt ddX_lt ddY_lt')
    X_rheel, Y_rheel, dX_rheel, dY_rheel, ddX_rheel, ddY_rheel = sp.symbols('X_rheel Y_rheel dX_rheel dY_rheel ddX_rheel ddY_rheel')

    # Defining torques
    tau_lt, tau_la, tau_lk, tau_lh, tau_hat, tau_rt, tau_ra, tau_rk, tau_rh = sp.symbols('tau_lt tau_la tau_lk tau_lh tau_hat tau_rt tau_ra tau_rk tau_rh')

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
    k_rt, k_ra, k_rk, k_rh = sp.symbols('k_rt k_ra k_rk k_rh')
    k_hat = sp.symbols('k_hat')
    C_lt, C_la, C_lk, C_lh = sp.symbols('C_lt C_la C_lk C_lh')
    C_rt, C_ra, C_rk, C_rh = sp.symbols('C_rt C_ra C_rk C_rh')
    C_hat = sp.symbols('C_hat')
    B = sp.symbols('B')

    # Parameters for Jacobian
    rs_f, rs_s, rs_t, rs_h = sp.symbols('rs_f rs_s rs_t rs_h')

    # Additional parameters
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

    # Left foot calculations
    VA__D_R_H = sp.Matrix([dX_lt, dY_lt, 0])
    aA = sp.Matrix([ddX_lt, ddY_lt, 0])
    Omega_lf = sp.Matrix([0, 0, dq1])
    alpha_lf = sp.Matrix([0, 0, ddq1])

    # Relative positions
    rG_lf = sp.Matrix([-r_f * sp.sin(q1), r_f * sp.cos(q1), 0])
    rBA = sp.Matrix([-L_f * sp.sin(q1), L_f * sp.cos(q1), 0])
    rIA = sp.Transpose(sp.Matrix([-rs_f * sp.sin(q1), rs_f * sp.cos(q1), 0]))

    # Velocity and acceleration
    VG_lf = VA__D_R_H + Omega_lf.cross(rG_lf)
    VB = VA__D_R_H + Omega_lf.cross(rBA)
    aG_lf = aA + alpha_lf.cross(rG_lf) + Omega_lf.cross(Omega_lf.cross(rG_lf))
    aB = aA + alpha_lf.cross(rBA) + Omega_lf.cross(Omega_lf.cross(rBA))

    XA__D_R_H = X_lt
    YA__D_R_H = Y_lt
    rA__D_R_H = sp.Transpose(sp.Matrix([XA__D_R_H, YA__D_R_H, 0]))

    X_lheel__D_R_H = XA__D_R_H - L_heelTotoe * sp.sin(sp.pi - alpha + q1)
    Y_lheel__D_R_H = YA__D_R_H + L_heelTotoe * sp.cos(sp.pi - alpha + q1)
    V_lheel__D_R_H = VA__D_R_H + Omega_lf.cross(sp.Matrix([-L_heelTotoe * sp.sin(sp.pi - alpha + q1), L_heelTotoe * sp.cos(sp.pi - alpha + q1), 0]))

    # Energy calculations
    T1 = 1/2 * m_f * (VG_lf.dot(VG_lf)) + 1/2 * I_f * (Omega_lf.dot(Omega_lf))
    V1 = m_f * g * (Y_lt + r_f * sp.cos(q1))
    Ve_lt = 1/2 * k_lt * (alpha - q1 - sp.pi/2)**2
    R_lt = 1/2 * C_lt * (Omega_lf.dot(Omega_lf))


    # GRF for Left toe
    # Penetration depth and rate (compression positive)
    delta_left_toe = sp.Max(-Y_lt, 0)
    contact_active = (1 + sp.tanh(delta_left_toe / eps_phi)) / 2
    vel_smooth = sp.tanh((-dY_lt) / v_eps)
    delta_dot_left_toe = (-dY_lt) * contact_active * vel_smooth
    # Hunt-Crossley normal force
    # F_n = k * delta^n * (1 + c * delta_dot)
    # The damping term is proportional to penetration depth (key feature!)
    Fy_left_toe_hc = k_hc * delta_left_toe**n_hc * (1 + c_hc * delta_dot_left_toe)
    # Stribeck friction curve for smoother transition
    mu_eff_lt = mu_d + (mu_s - mu_d) * sp.exp(-(dX_lt/v_s)**2)
    # Friction force with smooth saturation
    Fx_left_toe_hc = -mu_eff_lt * Fy_left_toe_hc * sp.tanh(dX_lt / eps_v)
    # Combined ground reaction force vector
    F_left_toe__D_R_H = sp.Matrix([Fx_left_toe_hc, Fy_left_toe_hc, 0])


    # Left shank
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
    V2 = m_s * g * (Y_lt + L_f * sp.cos(q1) + r_s * sp.cos(q2))
    Ve_la = 1/2 * k_la * (q2 - q1)**2
    R_la = 1/2 * C_la * ((Omega_ls - Omega_lf).dot(Omega_ls - Omega_lf))

    # Left thigh
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
    V3 = m_t * g * (Y_lt + L_f * sp.cos(q1) + L_s * sp.cos(q2) + r_t * sp.cos(q3))
    Ve_lk = 1/2 * k_lk * (q3 - q2)**2
    R_lk = 1/2 * C_lk * ((Omega_lt - Omega_ls).dot(Omega_lt - Omega_ls))

    # HAT (Head, Arms, Torso)
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
    V4 = m_h * g * (Y_lt + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) + r_h * sp.cos(q4))
    Ve_lh = 1/2 * k_lh * (q4 - q3)**2
    R_lh = 1/2 * C_lh * ((Omega_h - Omega_lt).dot(Omega_h - Omega_lt))

    # Right thigh
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
    V5 = m_t * g * (Y_lt + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) - s_t * sp.cos(q5))
    Ve_rh = 1/2 * k_rh * (q5 - q4)**2
    R_rh = 1/2 * C_rh * ((Omega_rt - Omega_h).dot(Omega_rt - Omega_h))

    # Right shank
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
    V6 = m_s * g * (Y_lt + L_f * sp.cos(q1) + L_s * sp.cos(q2) + L_t * sp.cos(q3) - L_t * sp.cos(q5) - s_s * sp.cos(q6))
    Ve_rk = 1/2 * k_rk * (q6 - q5)**2
    R_rk = 1/2 * C_rk * ((Omega_rs - Omega_rt).dot(Omega_rs - Omega_rt))

    # Positions of points B, C, D, E, F, G (left foot to right ankle)
    XB__D_R_H = XA__D_R_H - L_f * sp.sin(q1)
    YB__D_R_H = YA__D_R_H + L_f * sp.cos(q1)
    rB__D_R_H = sp.Transpose(sp.Matrix([XB__D_R_H, YB__D_R_H, 0]))

    XC__D_R_H = XB__D_R_H - L_s * sp.sin(q2)
    YC__D_R_H = YB__D_R_H + L_s * sp.cos(q2)
    rC__D_R_H = sp.Transpose(sp.Matrix([XC__D_R_H, YC__D_R_H, 0]))

    XD__D_R_H = XC__D_R_H - L_t * sp.sin(q3)
    YD__D_R_H = YC__D_R_H + L_t * sp.cos(q3)
    rD__D_R_H = sp.Transpose(sp.Matrix([XD__D_R_H, YD__D_R_H, 0]))

    XE__D_R_H = XD__D_R_H - L_h * sp.sin(q4)
    YE__D_R_H = YD__D_R_H + L_h * sp.cos(q4)
    rE__D_R_H = sp.Transpose(sp.Matrix([XE__D_R_H, YE__D_R_H, 0]))

    XF__D_R_H = XD__D_R_H + L_t * sp.sin(q5)
    YF__D_R_H = YD__D_R_H - L_t * sp.cos(q5)
    rF__D_R_H = sp.Transpose(sp.Matrix([XF__D_R_H, YF__D_R_H, 0]))

    XG__D_R_H = XF__D_R_H + L_s * sp.sin(q6)
    YG__D_R_H = YF__D_R_H - L_s * sp.cos(q6)
    rG__D_R_H = sp.Transpose(sp.Matrix([XG__D_R_H, YG__D_R_H, 0]))

    # Right foot
    # There is a cut in right ankle joint to simulate colsed chain with Lagrange method (see the article: https://doi.org/10.1016/j.mechmachtheory.2021.104294)
    # This cut makes the right foot seperate from other bodies
    # here the heel of right foot is always on the ground and new constraints are: X_rheel__D_R_H = Y_rheel__D_R_H = const (2l refers to "2" leg stance "left" leg at back)
    # q7 is equal to the q7 defined in left leg stance
    # Based on cutting joint method calculation is done for the right foot

    X_rheel, Y_rheel, dX_rheel, dY_rheel, ddX_rheel, ddY_rheel = sp.symbols('X_rheel Y_rheel dX_rheel dY_rheel ddX_rheel ddY_rheel')

    s_f = L_f - r_f
    ss_f = L_f - rs_f
    V_rheel__D_R_H = sp.Matrix([dX_rheel, dY_rheel, 0])
    a_rheel__D_R_H = sp.Matrix([ddX_rheel, ddY_rheel, 0])
    Omega_rf = sp.Matrix([0, 0, dq7])
    alpha_rf = sp.Matrix([0, 0, ddq7])

    # Positions of points G, H, and right heel calculated from cutted right foot
    X_rheel__D_R_H = X_rheel
    Y_rheel__D_R_H = Y_rheel

    XG_cutted = X_rheel__D_R_H + B * sp.cos(beta + sp.pi/2 - alpha + q7)
    YG_cutted = Y_rheel__D_R_H + B * sp.sin(beta + sp.pi/2 - alpha + q7)
    rG_cutted = sp.Matrix([XG_cutted, YG_cutted, 0])

    XH__D_R_H = XG_cutted + L_f * sp.sin(q7)
    YH__D_R_H = YG_cutted - L_f * sp.cos(q7)
    rH__D_R_H = sp.Matrix([XH__D_R_H, YH__D_R_H, 0])

    X_rbone__D_R_H = XH__D_R_H
    Y_rbone__D_R_H = YH__D_R_H

    # GRF for right heel
    # Penetration depth and rate (compression positive)
    delta_right_heel = sp.Max(-Y_rheel, 0)
    contact_active = (1 + sp.tanh(delta_right_heel / eps_phi)) / 2
    vel_smooth = sp.tanh((-dY_rheel) / v_eps)
    delta_dot_right_heel = (-dY_rheel) * contact_active * vel_smooth
    # Hunt-Crossley normal force
    # F_n = k * delta^n * (1 + c * delta_dot)
    # The damping term is proportional to penetration depth (key feature!)
    Fy_right_heel_hc = k_hc * delta_right_heel**n_hc * (1 + c_hc * delta_dot_right_heel)
    # Stribeck friction curve for smoother transition
    mu_eff_rh = mu_d + (mu_s - mu_d) * sp.exp(-(dX_rheel/v_s)**2)
    # Friction force with smooth saturation
    Fx_right_heel_hc = -mu_eff_rh * Fy_right_heel_hc * sp.tanh(dX_rheel / eps_v)
    # Combined ground reaction force vector
    F_right_heel__D_R_H = sp.Matrix([Fx_right_heel_hc, Fy_right_heel_hc, 0])


    # Relative positions for right foot
    rG_rf = sp.Matrix([X_rheel + B * sp.cos(beta + sp.pi/2 - alpha + q7) + s_f * sp.cos(sp.pi/2 - q7),
                       Y_rheel__D_R_H + B * sp.sin(beta + sp.pi/2 - alpha + q7) - s_f * sp.sin(sp.pi/2 - q7), 0])
    rHG = sp.Matrix([XG_cutted - XH__D_R_H, YG_cutted - YH__D_R_H, 0])
    rOG = sp.Transpose(sp.Matrix([ss_f * sp.cos(q7 - sp.pi/2), -ss_f * sp.sin(q7 - sp.pi/2), 0]))

    # Velocity and acceleration for right foot
    VG_cutted = V_rheel__D_R_H + Omega_rf.cross(rG_cutted)
    VG_rf = V_rheel__D_R_H + Omega_rf.cross(rG_rf)
    VH__D_R_H = V_rheel__D_R_H + Omega_rf.cross(rH__D_R_H)

    aG_cutted = a_rheel__D_R_H + alpha_rf.cross(rG_cutted) + Omega_rf.cross(Omega_rf.cross(rG_cutted))
    aG_rf = a_rheel__D_R_H + alpha_rf.cross(rG_rf) + Omega_rf.cross(Omega_rf.cross(rG_rf))
    aH = a_rheel__D_R_H + alpha_rf.cross(rH__D_R_H) + Omega_rf.cross(Omega_rf.cross(rH__D_R_H))

    # Energy calculations for right foot
    T7 = 1/2 * m_f * (VG_rf.dot(VG_rf)) + 1/2 * I_f * (Omega_rf.dot(Omega_rf))
    V7 = m_f * g * (YG_cutted - s_f * sp.sin(sp.pi/2 - q7))
    Ve_ra = 1/2 * k_ra * (q7 - q6)**2
    R_ra = 1/2 * C_ra * ((Omega_rf - Omega_rs).dot(Omega_rf - Omega_rs))

    # Ankle Constraint
    # In two leg stance we should introduce two new constraints that here are X and Y position of right heel
    # But in cutting joint method we change this constraints to XG = XG_cutted and YG = YG_cutted
    # So the constraints are making G point to be at same point in each of two seperated system of bodies
    # After creating the constraints their second derivative is needed in the dynamics equation which are called acceleration constraints
    # These acceleration constraints can be written in a matrix form  PHI : D*qdd=gamma
    # We need D and gamma matrices in the equations so we define them here
    # D is a 2*7 matrix and Gamma is a 2*1 matrix

    # Define constraints
    Phi_1__D_R_H = XG__D_R_H - XG_cutted
    Phi_2__D_R_H = YG__D_R_H - YG_cutted

    # Generalized coordinates and foot positions
    q = [q1, q2, q3, q4, q5, q6, q7]
    dq = [dq1, dq2, dq3, dq4, dq5, dq6, dq7]  # Only joint angle velocities (7 elements)
    foot_pos = [X_lt, Y_lt, X_rheel, Y_rheel]
    dfoot_pos = [dX_lt, dY_lt, dX_rheel, dY_rheel]  # Foot position velocities (4 elements)

    # Compute Jacobian (D) - first time derivatives
    D = sp.Matrix([[sp.diff(Phi_1__D_R_H, var) for var in q + foot_pos],
                   [sp.diff(Phi_2__D_R_H, var) for var in q + foot_pos]])

    # Compute dD/dq (derivative of D with respect to q)
    D_dq = sp.Matrix([[sp.diff(D[i, j], q[k]) for j in range(len(q + foot_pos))] for i in range(2) for k in range(len(q))])
    # D_dq is now 2x11x7, but we handle it row-wise for each constraint

    # Compute dD/dfoot_pos (derivative of D with respect to foot_pos)
    D_dfoot = sp.Matrix([[sp.diff(D[i, j], foot_pos[k]) for j in range(len(q + foot_pos))] for i in range(2) for k in range(len(foot_pos))])
    # D_dfoot is 2x11x4

    # Compute Gamma (2x1 matrix)
    Gamma = sp.Matrix.zeros(2, 1)
    for i in range(2):  # Two constraint equations
        for j in range(len(q)):  # First sum over joint angles q
            for k in range(len(q)):
                Gamma[i] -= sp.diff(D[i, j], q[k]) * dq[k] * dq[j]

        for j in range(len(q), len(q + foot_pos)):  # Second sum over foot positions
            foot_idx = j - len(q)  # Map to foot position index
            for k in range(len(foot_pos)):
                Gamma[i] -= sp.diff(D[i, j], foot_pos[k]) * dfoot_pos[k] * dfoot_pos[foot_idx]


    # Define symbols for generalized coordinates and their derivatives
    q = list(symbols('q1 q2 q3 q4 q5 q6 q7 X_lt Y_lt X_rheel Y_rheel'))
    dq = list(symbols('dq1 dq2 dq3 dq4 dq5 dq6 dq7 dX_lt dY_lt dX_rheel dY_rheel'))
    ddq = list(symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 ddX_lt ddY_lt ddX_rheel ddY_rheel'))

    # Lagrange’s Equations of Motion
    T = T1 + T2 + T3 + T4 + T5 + T6 + T7
    Ve = Ve_lt + Ve_la + Ve_lk + Ve_lh + Ve_rh + Ve_rk + Ve_ra
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
    dW = tau_lt * Omega_lf[2] + tau_la * (Omega_ls[2] - Omega_lf[2]) + tau_lk * (Omega_lt[2] - Omega_ls[2]) + \
         tau_lh * (Omega_h[2] - Omega_lt[2]) + tau_rh * (Omega_rt[2] - Omega_h[2]) + tau_rk * (Omega_rs[2] - Omega_rt[2]) + \
         tau_ra * (Omega_rf[2] - Omega_rs[2])

    # Rayleigh dissipation function
    R = R_lt + R_la + R_lk + R_lh + R_rh + R_rk + R_ra

    # Generalized forces Q
    Q = sp.Matrix([sp.diff(dW, dq_i) for dq_i in dq])

    # Jacobian of left toe and rgiht heel position with respect to generalized coordinates
    J_left_toe = sp.Matrix([XA__D_R_H, YA__D_R_H, 0]).jacobian(sp.Matrix(q))
    J_right_heel = sp.Matrix([X_rheel__D_R_H, Y_rheel__D_R_H, 0]).jacobian(sp.Matrix(q))

    # Total generalized ground force: toe + heel contributions
    Q_ground = J_left_toe.T * F_left_toe__D_R_H + J_right_heel.T * F_right_heel__D_R_H

    # Velocity-dependent generalized forces Q_v
    Q_v = -sp.Matrix([sp.diff(R, dq_i) for dq_i in dq])

    # Damping matrix C_v
    C_v = sp.Matrix([[sp.diff(Q_v_i, dq_j) for dq_j in dq] for Q_v_i in Q_v])

    # Total energy E
    E__D_R_H = T + V

    # Total generalized forces
    Q_total__D_R_H = Q + Q_ground

    # Time derivative of kinetic energy dE_dt
    dT_dt_ls = sum(sp.Matrix([sp.diff(T, q_i) * dq_i + sp.diff(T, dq_i) * ddq_i for q_i, dq_i, ddq_i in zip(q, dq, ddq)]))

    # Time derivative of energy dE_dt
    dV_dt_ls = sum(sp.Matrix([sp.diff(V, q_i) * dq_i + sp.diff(V, dq_i) * ddq_i for q_i, dq_i, ddq_i in zip(q, dq, ddq)]))

    # Time derivative of energy dE_dt
    dE_dt__D_R_H = sum(sp.Matrix([sp.diff(E__D_R_H, q_i) * dq_i + sp.diff(E__D_R_H, dq_i) * ddq_i for q_i, dq_i, ddq_i in zip(q, dq, ddq)]))

    # Work done by non-conservative forces dW_nc
    dW_nc = sum([sp.diff(Q_i, dq_j) * dq_j for Q_i, dq_j in zip(Q, dq)]) - sum([sp.diff(R, dq_i) * dq_i for dq_i in dq])


    return {
    "name": "D_R_H",

    # Generalized coordinates
    "q": q,
    "dq": dq,
    "ddq": ddq,

    # Dynamics
    "M": M,
    "C": C + C_v,
    "G": G,
    "D": D,
    "Gamma": Gamma,
    "F_generalized": Q_ground,
    "Phi_1": Phi_1__D_R_H,
    "Phi_2": Phi_2__D_R_H,

    # Contact forces (raw)
    "contact": {
        "left_toe_force": F_left_toe__D_R_H,
        "right_heel_force": F_right_heel__D_R_H
    },

    # Positions (symbolic, BEFORE substitution)
    "positions": {
        "A": (XA__D_R_H, YA__D_R_H),
        "B": (XB__D_R_H, YB__D_R_H),
        "C": (XC__D_R_H, YC__D_R_H),
        "D": (XD__D_R_H, YD__D_R_H),
        "E": (XE__D_R_H, YE__D_R_H),
        "F": (XF__D_R_H, YF__D_R_H),
        "G": (XG__D_R_H, YG__D_R_H),
        "H": (XH__D_R_H, YH__D_R_H),
        "lheel": (X_lheel__D_R_H, Y_lheel__D_R_H),
        "rheel": (X_rheel__D_R_H, Y_rheel__D_R_H)
    },

    # Velocities
    "velocities": {
        "A": VA__D_R_H,
        "lheel": V_lheel__D_R_H,
        "H": VH__D_R_H,
        "rheel": V_rheel__D_R_H
    }
}
