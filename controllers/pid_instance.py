from controllers.pid import PIDController


def create_pid_controllers():
    """
    Create PID controllers for single and double support.
    Phase is NOT stored inside the controller anymore.
    """

    return PIDController(
           Kp=[400, 700, 700, 3300, 700, 500, 400],
           Kd=[20, 40, 50, 70, 50, 40, 20],
           Ki=[0.1] * 7)
