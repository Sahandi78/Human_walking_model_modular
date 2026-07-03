import sympy as sp
import math
from sympy import cos

# Define symbolic variables (ensure these match those used in your equations)
g, h, m_to, m_f, m_s, m_t, m_h, L_f, L_s, L_t, L_h, L_to, r_f, r_s, r_t, r_h, r_to = sp.symbols('g h m_to m_f m_s m_t m_h L_f L_s L_t L_h L_to r_f r_s r_t r_h r_to')
I_f, I_s, I_t, I_h, I_to, B, alpha, beta, L_heelTotoe, M_human_tot = sp.symbols('I_f I_s I_t I_h I_to B alpha beta L_heelTotoe M_human_tot')
k_lt, k_rt, k_la, k_ra, k_lk, k_rk, k_lh, k_rh, C_lt, C_rt, C_la, C_ra, C_lk, C_rk, C_lh, C_rh = sp.symbols('k_lt k_rt k_la k_ra k_lk k_rk k_lh k_rh C_lt C_rt C_la C_ra C_lk C_rk C_lh C_rh')
k_lbone, k_lheel, k_rbone, k_rheel = sp.symbols('k_lbone k_lheel k_rbone k_rheel')
C_lbone, C_lheel, C_rbone, C_rheel = sp.symbols('C_lbone C_lheel C_rbone C_rheel')
n_hc, k_hc, c_hc = sp.symbols('n_hc k_hc c_hc')
mu_s, mu_d, v_s = sp.symbols('mu_s mu_d v_s')
eps_v = sp.symbols('eps_v')
eps_phi = sp.symbols('eps_phi')
v_eps = sp.symbols('v_eps')



values_to_substitute = {
g : 9.81,
m_to : 0.15, # Toe mass
m_f : 0.0145*56.7, # Foot mass
m_s : 0.0465*56.7, # Shank mass
m_t : 0.1*56.7, # Thigh mass
m_h : 0.678*56.7, # HAT mass
L_s : 0.425, # Length of shank
L_t : 0.314, # Length of thigh
L_h : 0.709152, # Length of HAT
L_to : 0.051, # Length of toe
L_f : 0.121, # Length of foot
r_f : 0.121*0.5, # COM length of foot
r_s : 0.433*0.425, # COM length of shank
r_t : 0.433*0.314, # COM length of thigh
r_h : 0.626*0.709152, # COM length of HAT
r_to : 0.051/2, # COM length of toe
I_f : (0.475*0.121)**2 * 0.0145*56.7, # Inertia of foot
I_s : (0.302*0.425)**2 * 0.0465*56.7, # Inertia of shank
I_t : (0.323*0.314)**2 * 0.1*56.7, # Inertia of thihg
I_h : (0.496*0.709152)**2 * 0.678*56.7, # # Inertia of HAT
I_to : 0.005, # Inertia of toe
B : 0.075, # heel to ankle line length
alpha : 2.621, # toe to ankle line angle 
beta : 0.967, # heel to ankle line angle
L_heelTotoe : 0.121*cos(math.pi-2.621)+0.075*cos(0.967), # Heel to toe length
h: 0.0001, # Constant for contact detection (in meter)
k_hc: 5e5,      # Hunt-crossley stiffness (N/m)
c_hc: 2e-1,      # Hunt-crossley damping (N.s/m)
n_hc: 1.5,      # Standard Hertzian 
eps_v: 1e-2,   
eps_phi: 1e-4, 
v_eps: 0.01,
mu_s: 0.9,      # static friction
mu_d: 0.7,      # dynamic friction
v_s: 0.01,       # Stribeck velocity
k_lt: 0,
k_rt: 0,
k_ra: 0,
k_la: 0,
k_rk: 0,
k_lk: 0,
k_rh: 0,
k_lh: 0,
C_lt: 0,
C_rt: 0,
C_ra: 0,
C_la: 0,
C_rk: 0,
C_lk: 0,
C_rh: 0,
C_lh: 0
}
