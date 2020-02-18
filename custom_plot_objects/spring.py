import matplotlib.pyplot as plt
import numpy as np

class Spring:
    """Plot a spring"""
    def __init__(self,p1,p2,ax = None):
        if ax is None:
            ax = plt.gca()
        self.springObj, = ax.plot([], [], 'k-',lw=2)

        x1,y1 = p1
        x2,y2 = p2

        L = np.sqrt((x2-x1)**2+(y2-y1)**2)
        theta = np.arctan2((y2-y1),(x2-x1))

        self.ReplotSpring(p1,theta,L)

    @classmethod
    def FromLengthAndAngle(cls, theta,L,p1=(0,0),ax = None):
        x1,y1 = p1
        x2 = x1+np.cos(theta)*L
        y2 = y1 + np.sin(theta)*L
        return cls(p1,(x2,y2),ax)
    


    def ReplotSpring(self,p1,theta,L):
        x1,y1 = p1
        # Spring turn radius, number of turns
        rs, ns = 0.05, 10
        # Number of data points for the helix
        Ns = 1000
        # We don't draw coils all the way to the end of the pendulum:
        # pad a bit from the anchor and from the bob by these number of points
        ipad1, ipad2 = 100, 150

        # calculate spring length and angle
        w = np.linspace(0, L, Ns)

        xp = np.zeros(Ns)
        xp[ipad1:-ipad2] = rs * np.sin(2*np.pi * ns * w[ipad1:-ipad2] / L)
        # ... then rotate it to align with  the pendulum and plot.
        R = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
        xs, ys = - R @ np.vstack((xp, w))
        self.springObj.set_data(xs+x1, ys+y1)



