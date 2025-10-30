import numpy as np
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def add_circular_colorbar(ax, mappable, orientation="vertical", fraction=0.15, pad=0.1):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax,
                       width="5%", 
                       height="50%",
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0.2, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    cbar = plt.colorbar(mappable, cax=axins, orientation=orientation, ticks=[0, 180])
    cbar.ax.set_yticklabels(['0°', '180°'])  
    return cbar

# Create a black canvas
canvas = np.zeros((512, 512), dtype=np.uint8)

# Draw random ellipses
for _ in range(10):
    center = (np.random.randint(50, 462), np.random.randint(50, 462))
    axes = (np.random.randint(20, 50), np.random.randint(10, 20))
    angle = np.random.randint(0, 180)
    startAngle = 0
    endAngle = 360
    color = 255 # white
    thickness = -1 # fill the ellipse

    cv2.ellipse(canvas, center, axes, angle, startAngle, endAngle, color, thickness)

# Detect ellipses from the canvas
contours, _ = cv2.findContours(canvas, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

fig, ax = plt.subplots()
ax.imshow(canvas, cmap='gray')

colormap = plt.cm.hsv

for cnt in contours:
    ellipse = cv2.fitEllipse(cnt)
    (x0, y0), (minor_axis, major_axis), angle = ellipse
    
    # Map the orientation angle to a colormap
    color = colormap(angle/180.)  # normalize angle to [0,1] to match colormap

    ellipse_patch = mpatches.Ellipse(
        (x0, y0),
        major_axis,
        minor_axis,
        angle=angle,
        fill=False,
        edgecolor=color,
        linewidth=2
    )
    ax.add_patch(ellipse_patch)

# Create a fake image for the colormap
img = np.array([[0, 180]])
mappable = ax.imshow(img, cmap=colormap, aspect='auto', visible=False)
add_circular_colorbar(ax, mappable, orientation="vertical")

plt.show()

