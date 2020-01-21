import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def xyline(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

x = np.linspace(0, 3 * np.pi, 500)
y = np.sin(x)
dydx = np.cos(0.5 * (x[:-1] + x[1:]))  # first derivative

g = np.linspace(0.1, 0.7, 10)
colors = np.array([0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4])
rs = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0])

# Create a set of line segments so that we can color them individually
# This creates the points as a N x 1 x 2 array so that we can stack points
# together easily to get the segments. The segments array for line collection
# needs to be (numlines) x (points per line) x 2 (for x and y)
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
points = np.array([rs, g]).T.reshape(-1, 1, 2)
s2 = np.concatenate([points[:-1], points[1:]], axis=1)


fig = plt.figure()
ax = fig.add_subplot(111)

norm = plt.Normalize(dydx.min(), dydx.max())

lc = LineCollection(segments, cmap='Reds', norm=norm)
lc.set_array(dydx)
lc.set_linewidth(2)
line = ax.add_collection(lc)

lc = LineCollection(s2, cmap='Reds', norm=norm)
lc.set_array(colors)
lc.set_linewidth(2)
ll2 = ax.add_collection(lc)
fig.colorbar(ll2, ax=ax)

ax.set_xlim(-0.2,6)
ax.set_ylim(g.min() - 1, g.max() + 1)

'''
# Create a continuous norm to map from data points to colors
norm = plt.Normalize(dydx.min(), dydx.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
# Set the values used for colormapping
lc.set_array(dydx)
lc.set_linewidth(2)
line = axs[0].add_collection(lc)
fig.colorbar(line, ax=axs[0])

# Use a boundary norm instead
cmap = ListedColormap(['r', 'g', 'b'])
norm = BoundaryNorm([-1, -0.5, 0.5, 1], cmap.N)
lc = LineCollection(segments, cmap=cmap, norm=norm)
lc.set_array(dydx)
lc.set_linewidth(2)
line = axs[1].add_collection(lc)
fig.colorbar(line, ax=axs[1])
'''
