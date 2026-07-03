# plotting/animation.py

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


def animate_walking(results, step=50):

    # Extract data
    XA = results["XA"]; YA = results["YA"]
    XB = results["XB"]; YB = results["YB"]
    XC = results["XC"]; YC = results["YC"]
    XD = results["XD"]; YD = results["YD"]
    XE = results["XE"]; YE = results["YE"]
    XF = results["XF"]; YF = results["YF"]
    XG = results["XG"]; YG = results["YG"]
    XH = results["XH"]; YH = results["YH"]

    X_lheel = results["X_lheel"]
    Y_lheel = results["Y_lheel"]
    X_rheel = results["X_rheel"]
    Y_rheel = results["Y_rheel"]

    fig, ax = plt.subplots(figsize=(8, 6))

    # -------------------------------------------------
    # Init function (important for blitting stability)
    # -------------------------------------------------
    def init():
        ax.set_xlim(-1, 2)
        ax.set_ylim(-0.5, 2)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.axis('equal')
        return []

    # -------------------------------------------------
    # Update function
    # -------------------------------------------------
    def update(i):

        ax.clear()

        # Ground
        ax.axhline(0, color='black', linestyle='--', linewidth=1)

        # Limb segments
        ax.plot([XB[i], XC[i]], [YB[i], YC[i]], 'r-', linewidth=2)
        ax.plot([XC[i], XD[i]], [YC[i], YD[i]], 'g-', linewidth=2)
        ax.plot([XD[i], XE[i]], [YD[i], YE[i]], 'b-', linewidth=5)
        ax.plot([XD[i], XF[i]], [YD[i], YF[i]], 'g-', linewidth=2)
        ax.plot([XF[i], XG[i]], [YF[i], YG[i]], 'r-', linewidth=2)

        # Feet (triangles)
        Xl = [XA[i], XB[i], X_lheel[i]]
        Yl = [YA[i], YB[i], Y_lheel[i]]

        Xr = [XG[i], XH[i], X_rheel[i]]
        Yr = [YG[i], YH[i], Y_rheel[i]]

        ax.fill(Xl, Yl, 'r')
        ax.fill(Xr, Yr, 'r')

        # Joint markers
        ax.plot(
            [XB[i], XC[i], XD[i], XF[i], XG[i]],
            [YB[i], YC[i], YD[i], YF[i], YG[i]],
            'go'
        )

        ax.set_xlim(-1, 2)
        ax.set_ylim(-0.5, 2)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title("Walking Animation")
        ax.axis('equal')

        return []

    ani = FuncAnimation(
        fig,
        update,
        frames=np.arange(0, len(XA), step),
        init_func=init,
        blit=True,
        repeat=True
    )

    plt.show()
