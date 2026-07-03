# plotting/torque_plots.py

import numpy as np
import matplotlib.pyplot as plt


def plot_joint_torques(results):
    """
    Plot joint torques using provided compute_torques function.
    """

    torques = results["Joint_torques"]
    time = results["time"]

    plt.figure(figsize=(12, 8))
    for i in range(7):
        plt.subplot(3, 3, i + 1)
        plt.plot(time, torques[:, i])
        plt.title(f'Joint {i+1}')
        plt.grid(True)

    plt.tight_layout()
    plt.show()
