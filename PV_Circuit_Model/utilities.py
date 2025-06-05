import numpy as np
import numbers
from matplotlib import pyplot as plt
import matplotlib.patches as patches

def interp_(x, xp, fp):
    """
    Like np.interp, but extrapolates linearly outside the bounds using the slopes
    of the first two and last two points.
    """
    xp_ = xp.copy()
    fp_ = fp.copy()
    is_number = False
    if isinstance(x, numbers.Number):
        x = np.array([x])
        is_number = True
    x_ = x.copy()
    while xp_[0]==xp_[1]:
        xp_ = xp_[1:]
        fp_ = fp_[1:]
    while xp_[-1]==xp_[-2]:
        xp_ = xp_[:-1]
        fp_ = fp_[:-1]
    if xp_[0] > xp_[-1]:
        xp_ = -1*xp_
        x_ = -1*x_
    y_multiplier = 1.0
    if fp_[0] > fp_[-1]:
        y_multiplier = -1.0
        fp_ = -1*fp_
    slope_left = (fp_[1] - fp_[0]) / (xp_[1] - xp_[0])
    slope_right = (fp_[-1] - fp_[-2]) / (xp_[-1] - xp_[-2])
    y = np.interp(x_, xp_, fp_)
    y[x_ < xp_[0]] = fp_[0] + slope_left * (x_[x_ < xp_[0]] - xp_[0])
    y[x_ > xp_[-1]] = fp_[-1] + slope_right * (x_[x_ > xp_[-1]] - xp_[-1])
    y *= y_multiplier
    if is_number:
        return y[0]
    return y

def draw_symbol(draw_func, ax=None,  x=0, y=0, color="black", text=None, **kwargs):
    draw_immediately = False
    if ax is None:
        draw_immediately = True
        _, ax = plt.subplots()
    draw_func(ax=ax, x=x, y=y, color=color, **kwargs)
    if text is not None:
        text_x = 0.14
        text_y = 0.0
        if draw_func==draw_CC_symbol:
            text_x = 0.21
        elif draw_func==draw_resistor_symbol:
            text_y = -0.15
        ax.text(x+text_x,y+text_y,text, va='center', fontsize=6, color=color)
    if draw_immediately:
        plt.show()

def draw_diode_symbol(ax=None, x=0, y=0, color="black", up_or_down="down", is_LED=False):
    dir = 1
    if up_or_down == "up":
        dir = -1
    if is_LED:
        circle = patches.Circle((x,y), 0.17, edgecolor=color,facecolor='white',linewidth=1.5, fill=False)
        ax.add_patch(circle)
    ax.arrow(x, y+0.075*dir, 0, -0.001*dir, head_width=0.15, head_length=0.15, color=color, fc=color, ec=color)
    line = plt.Line2D([x-0.075,x+0.075], [y-0.08*dir,y-0.08*dir], color=color, linewidth=2)
    ax.add_line(line)
    if is_LED:
        ax.arrow(x-0.05, y-0.05*dir, -0.15, -0.15*dir, head_width=0.05, head_length=0.05, fc='orange', ec='orange')
        ax.arrow(x-0.075, y+0.025*dir, -0.15, -0.15*dir, head_width=0.05, head_length=0.05, fc='orange', ec='orange')
    line = plt.Line2D([x,x], [y+0.08, y+0.4], color="black", linewidth=1.5)
    ax.add_line(line)
    line = plt.Line2D([x,x], [y-0.08, y-0.4], color="black", linewidth=1.5)
    ax.add_line(line)

def draw_forward_diode_symbol(ax, x=0, y=0, color="black"):
    draw_diode_symbol(ax=ax, x=x, y=y, color=color, up_or_down="down", is_LED=False)

def draw_reverse_diode_symbol(ax, x=0, y=0, color="black"):
    draw_diode_symbol(ax=ax, x=x, y=y, color=color, up_or_down="up", is_LED=False)

def draw_LED_diode_symbol(ax, x=0, y=0, color="black"):
    draw_diode_symbol(ax=ax, x=x, y=y, color=color, up_or_down="down", is_LED=True)

def draw_CC_symbol(ax, x=0, y=0, color="black"):
    circle = patches.Circle((x, y), 0.17, edgecolor=color,facecolor="white",linewidth=2, fill=True)
    ax.add_patch(circle)
    ax.arrow(x, y-0.12, 0, 0.14, head_width=0.1, head_length=0.1, width=0.01,  fc=color, ec=color)
    line = plt.Line2D([x,x], [y+0.18, y+0.4], color="black", linewidth=1.5)
    ax.add_line(line)
    line = plt.Line2D([x,x], [y-0.18, y-0.4], color="black", linewidth=1.5)
    ax.add_line(line)

def draw_resistor_symbol(ax, x=0, y=0, color="black"):
    ax.arrow(x, y-0.25, 0, 0.12, head_width=0.1, head_length=0.1, width=0.01,  fc="white", ec="white")
    dx = 0.075
    dy = 0.02
    ystart = y + 0.15
    line = plt.Line2D([x,x], [y+0.15, y+0.4], color="black", linewidth=1.5)
    ax.add_line(line)
    line = plt.Line2D([x,x], [y-0.09, y-0.4], color="black", linewidth=1.5)
    ax.add_line(line)
    for _ in range(3):
        line = plt.Line2D([x,x+dx], [ystart,ystart-dy], color=color, linewidth=1.5)
        ax.add_line(line)
        line = plt.Line2D([x+dx,x-dx], [ystart-dy,ystart-3*dy], color=color, linewidth=1.5)
        ax.add_line(line)
        line = plt.Line2D([x-dx,x], [ystart-3*dy,ystart-4*dy], color=color, linewidth=1.5)
        ax.add_line(line)
        ystart -= 4*dy

def draw_earth_symbol(ax, x=0, y=0, color="black"):
    for i in range(3):
        line = plt.Line2D([x-0.03*(i+1),x+0.03*(i+1)], [y+0.05*i,y+0.05*i], linewidth=2, color=color)
        ax.add_line(line)

def draw_pos_terminal_symbol(ax, x=0, y=0, color="black"):
    circle = patches.Circle((x, y), 0.04, edgecolor=color,facecolor="white",linewidth=2, fill=True)
    ax.add_patch(circle)

class RandomNumberGenerator():
    pass

class CappedAbsGaussian(RandomNumberGenerator):
    def __init__(self,mean,stdev,cap=None):
        self.mean = mean
        self.stdev = stdev
        self.cap = cap
    def generate(self, sample_size = 1):
        x = np.abs(np.random.normal(self.mean,self.stdev,sample_size))
        if self.cap is not None:
            x = np.minimum(x,self.cap)
        return x
