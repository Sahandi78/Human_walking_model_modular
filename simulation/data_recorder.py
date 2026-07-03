import numpy as np

class DataRecorder:
    """
    Class for recording all simulation data
    """

    def __init__(self, numeric_phases):
        # Time
        self.time = []

        # Joint states
        self.positions = []
        self.velocities = []

        # Phases
        self.phases = []
        self.phases_one_leg = []
        self.phases_two_leg = []

        # Lambda values for double stance
        self.lambdas = []

        # Foot kinematics
        self.foot_pos = []
        self.dfoot_pos = []

        # Ground reaction forces (Fx, Fy, Fz)
        self.GRF_left_toe = []
        self.GRF_left_heel = []
        self.GRF_right_toe = []
        self.GRF_right_heel = []

        # Reference trajectories
        self.q_ref = []
        self.dq_ref = []
        self.ddq_ref = []

        # Time step record
        self.dt_all = []
        
        # lambdify function
        self.numeric_phases = numeric_phases
        
        # Position and velcities 
        self.XA = []
        self.YA = []
        self.VA = []
        self.XB = []
        self.YB = []
        self.XC = []
        self.YC = []
        self.XD = []
        self.YD = []
        self.XE = []
        self.YE = []
        self.XF = []
        self.YF = []
        self.XG = []
        self.YG = []
        self.XH = []
        self.YH = []
        self.VH = []
        self.X_rheel = []
        self.Y_rheel = []
        self.X_lheel = []
        self.Y_lheel = []
        self.V_lheel = []
        self.V_rheel = []
        
        # Joint torques
        self.Joint_torques = []


    # ---------------------------------------------------------
    # Utility: ensure GRF is stored as a clean (3,) array
    # ---------------------------------------------------------
    @staticmethod
    def _to_vector(val, dim=2):
        """
        Convert GRF output to a numpy vector of length = dim.
        Pads with zeros if needed.
        """
        if val is None:
            return np.zeros(dim)

        arr = np.asarray(val, dtype=float).squeeze()

        # Scalar → expand
        if arr.ndim == 0:
            arr = np.array([arr])

        # Pad if smaller
        if arr.shape[0] < dim:
            padded = np.zeros(dim)
            padded[:arr.shape[0]] = arr
            return padded

        # Truncate if larger
        if arr.shape[0] > dim:
            return arr[:dim]

        return arr
        
    
    # ---------------------------------------------------------
    # Utility: extract position and velocity of a point
    # ---------------------------------------------------------
    def get_point_state(self, label, t, y, current_phase):
        """
        Returns (X, Y, dX, dY) for a named point
        """
        numeric_phase = self.numeric_phases[current_phase]
        
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


    # ---------------------------------------------------------
    # Main record function
    # ---------------------------------------------------------
    def record(
        self,
        t,
        dt,
        y_current,
        current_phase,
        GRF,
        is_double_stance,
        torques,
        q_ref=None,
        dq_ref=None,
        ddq_ref=None
    ):
        """
        Record simulation data at a single time step.

        GRF = (F_lt, F_lh, F_rt, F_rh)
        Each force is a 3D vector [Fx, Fy, Fz]
        """

        F_lt, F_lh, F_rt, F_rh = GRF

        # -----------------------------
        # Time
        # -----------------------------
        self.time.append(t)
        self.dt_all.append(dt)

        # -----------------------------
        # Joint states
        # -----------------------------
        self.positions.append(y_current[:7])
        self.velocities.append(y_current[7:14])

        # -----------------------------
        # Phase & constraints
        # -----------------------------
        self.phases.append(current_phase)

        if is_double_stance:
            self.phases_two_leg.append(current_phase)
            self.lambdas.append(y_current[-2:])
            self.foot_pos.append(y_current[14:18])
            self.dfoot_pos.append(y_current[18:22])
        else:
            self.phases_one_leg.append(current_phase)
            self.foot_pos.append(y_current[14:16])
            self.dfoot_pos.append(y_current[16:18])

        # -----------------------------
        # Ground reaction forces (VECTOR)
        # -----------------------------
        self.GRF_left_toe.append(self._to_vector(F_lt))
        self.GRF_left_heel.append(self._to_vector(F_lh))
        self.GRF_right_toe.append(self._to_vector(F_rt))
        self.GRF_right_heel.append(self._to_vector(F_rh))

        # -----------------------------
        # Reference trajectories
        # -----------------------------
        if q_ref is not None:
            self.q_ref.append(q_ref)
        if dq_ref is not None:
            self.dq_ref.append(dq_ref)
        if ddq_ref is not None:
            self.ddq_ref.append(ddq_ref)
    
        # -----------------------------
        # Postions and velocities of different points
        # -----------------------------
        numeric_phase = self.numeric_phases[current_phase]
        
        XA_fun, YA_fun = numeric_phase["positions"]["A"]
        self.XA.append(XA_fun(t, y_current))
        self.YA.append(YA_fun(t, y_current))
        
        XB_fun, YB_fun = numeric_phase["positions"]["B"]
        self.XB.append(XB_fun(t, y_current))
        self.YB.append(YB_fun(t, y_current))
        
        XC_fun, YC_fun = numeric_phase["positions"]["C"]
        self.XC.append(XC_fun(t, y_current))
        self.YC.append(YC_fun(t, y_current))
        
        XD_fun, YD_fun = numeric_phase["positions"]["D"]
        self.XD.append(XD_fun(t, y_current))
        self.YD.append(YD_fun(t, y_current))
        
        XE_fun, YE_fun = numeric_phase["positions"]["E"]
        self.XE.append(XE_fun(t, y_current))
        self.YE.append(YE_fun(t, y_current))
        
        XF_fun, YF_fun = numeric_phase["positions"]["F"]
        self.XF.append(XF_fun(t, y_current))
        self.YF.append(YF_fun(t, y_current))
        
        XG_fun, YG_fun = numeric_phase["positions"]["G"]
        self.XG.append(XG_fun(t, y_current))
        self.YG.append(YG_fun(t, y_current))
        
        XH_fun, YH_fun = numeric_phase["positions"]["H"]
        self.XH.append(XH_fun(t, y_current))
        self.YH.append(YH_fun(t, y_current))

        X_lheel_fun, Y_lheel_fun = numeric_phase["positions"]["lheel"]
        self.X_lheel.append(X_lheel_fun(t, y_current))
        self.Y_lheel.append(Y_lheel_fun(t, y_current))
        
        X_rheel_fun, Y_rheel_fun = numeric_phase["positions"]["rheel"]
        self.X_rheel.append(X_rheel_fun(t, y_current))
        self.Y_rheel.append(Y_rheel_fun(t, y_current))
        
        self.VA.append(numeric_phase["velocities"]["A"](t, y_current))
        self.VH.append(numeric_phase["velocities"]["H"](t, y_current))
        self.V_lheel.append(numeric_phase["velocities"]["lheel"](t, y_current))
        self.V_rheel.append(numeric_phase["velocities"]["rheel"](t, y_current))
        
        # -----------------------------
        # Joint torques
        # -----------------------------
        self.Joint_torques.append(torques)

    # ---------------------------------------------------------
    # Export results
    # ---------------------------------------------------------
    def get_results(self):
        return {
            "time": np.array(self.time),
            "dt": np.array(self.dt_all),
            "positions": np.array(self.positions),
            "velocities": np.array(self.velocities),
            "phases": self.phases,
            "phases_one_leg": self.phases_one_leg,
            "phases_two_leg": self.phases_two_leg,
            #"lambdas": np.array(self.lambdas),
            #"foot_pos": np.array(self.foot_pos),
            #"dfoot_pos": np.array(self.dfoot_pos),

            # GRFs → (N, 3)
            "GRF_left_toe": np.vstack(self.GRF_left_toe),
            "GRF_left_heel": np.vstack(self.GRF_left_heel),
            "GRF_right_toe": np.vstack(self.GRF_right_toe),
            "GRF_right_heel": np.vstack(self.GRF_right_heel),

            "q_ref": np.array(self.q_ref) if self.q_ref else None,
            "dq_ref": np.array(self.dq_ref) if self.dq_ref else None,
            "ddq_ref": np.array(self.ddq_ref) if self.ddq_ref else None,
            
            "XA": np.array(self.XA),
            "YA": np.array(self.YA),
            "XB": np.array(self.XB),
            "YB": np.array(self.YB),
            "XC": np.array(self.XC),
            "YC": np.array(self.YC),
            "XD": np.array(self.XD),
            "YD": np.array(self.YD),
            "XE": np.array(self.XE),
            "YE": np.array(self.YE),
            "XF": np.array(self.XF),
            "YF": np.array(self.YF),
            "XG": np.array(self.XG),
            "YG": np.array(self.YG),
            "XH": np.array(self.XH),
            "YH": np.array(self.YH),
            "X_rheel": np.array(self.X_rheel),
            "Y_rheel": np.array(self.Y_rheel),
            "X_lheel": np.array(self.X_lheel),
            "Y_lheel": np.array(self.Y_lheel),
            "VA": np.array(self.VA),
            "VH": np.array(self.VH),
            "V_rheel": np.array(self.V_rheel),
            "V_lheel": np.array(self.V_lheel),
            
            "Joint_torques": np.array(self.Joint_torques),
        }