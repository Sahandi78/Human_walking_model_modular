# plotting/tracking_plots.py

import numpy as np
import matplotlib.pyplot as plt
import math


def plot_tracking(results, plot_q=7):
    """
    Plot joint position and velocity tracking.
    """

    positions = results["positions"]
    velocities = results["velocities"]
    q_ref = results["q_ref"]
    dq_ref = results["dq_ref"]

    if q_ref is None or dq_ref is None:
        print("Reference trajectories not available.")
        return

    ncols = 3
    nrows = math.ceil(plot_q / ncols)

    t = np.arange(len(positions))

    # ======================
    # Positions
    # ======================
    fig_pos, axes_pos = plt.subplots(nrows, ncols, figsize=(12, 8), sharex=True)
    axes_pos = axes_pos.flatten()
    fig_pos.suptitle('Joint Position Tracking', fontsize=16, fontweight='bold')

    for i in range(plot_q):
        ax = axes_pos[i]
        ax.plot(t, positions[:, i], label='Sim')
        ax.plot(t, q_ref[:, i], '--', label='Ref')
        ax.set_title(f'Joint {i+1}')
        ax.grid(True)
        ax.legend()

    for j in range(plot_q, len(axes_pos)):
        fig_pos.delaxes(axes_pos[j])

    plt.tight_layout()
    plt.show()

    # ======================
    # Velocities
    # ======================
    fig_vel, axes_vel = plt.subplots(nrows, ncols, figsize=(12, 8), sharex=True)
    axes_vel = axes_vel.flatten()
    fig_vel.suptitle('Joint Velocity Tracking', fontsize=16, fontweight='bold')

    for i in range(plot_q):
        ax = axes_vel[i]
        ax.plot(t, velocities[:, i], label='Sim')
        ax.plot(t, dq_ref[:, i], '--', label='Ref')
        ax.set_title(f'Joint {i+1}')
        ax.grid(True)
        ax.legend()

    for j in range(plot_q, len(axes_vel)):
        fig_vel.delaxes(axes_vel[j])

    plt.tight_layout()
    plt.show()
