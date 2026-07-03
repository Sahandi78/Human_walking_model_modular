"""
Numerical thresholds for gait event detection
"""

# Ground height (usually 0, but configurable)
GROUND_HEIGHT = 0.0

# Height tolerance (meters)
H_TOL = 0.001       # 1 mm

# Velocity thresholds (m/s)
V_STRIKE = -0.01    # downward
V_OFF = 0.01        # upward

# Contact force threshold (N)
FORCE_OFF = 5.0
