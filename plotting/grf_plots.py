# plotting/grf_plots.py

import numpy as np
import matplotlib.pyplot as plt


def plot_grf_2x2(results):
    """
    Publication-style 2x2 GRF plot.
    """

    time = results["time"]

    GRF_left_toe = results["GRF_left_toe"][:, 1]
    GRF_left_heel = results["GRF_left_heel"][:, 1]
    GRF_right_toe = results["GRF_right_toe"][:, 1]
    GRF_right_heel = results["GRF_right_heel"][:, 1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    axes[0, 0].plot(time, GRF_left_toe)
    axes[0, 0].set_title("Left Toe")

    axes[0, 1].plot(time, GRF_left_heel)
    axes[0, 1].set_title("Left Heel")

    axes[1, 0].plot(time, GRF_right_toe)
    axes[1, 0].set_title("Right Toe")

    axes[1, 1].plot(time, GRF_right_heel)
    axes[1, 1].set_title("Right Heel")

    for ax in axes.flat:
        ax.grid(True)

    plt.tight_layout()
    plt.show()
