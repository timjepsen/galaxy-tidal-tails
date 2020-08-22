#This program handles plotting an animation of the ODE results
#It's based on based on https://matplotlib.org/3.1.1/gallery/widgets/slider_demo.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button, RadioButtons
from Stage_2_Calculation import run
from Parameters import max_time

animation_running = False

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25)

plt.axis([-25, 25, -25, 25])

l_plot, = plt.plot([], [], ls = " ", marker = ".")
ax.margins(x=0)

axcolour = 'lightgoldenrodyellow'

odeSolution = run()

def prepare_to_plot_instantaneous_positions(timeIndex):
    """Takes the data from SolveOde() and gets it ready to be plotted (by visualising every particle at a specific point
    in time)"""
    x_coords = []
    y_coords = []
    z_coords = []

    for i in range(0, len(odeSolution.y), 6):
        x_coords.append(odeSolution.y[i][timeIndex])
        y_coords.append(odeSolution.y[i+1][timeIndex])
        z_coords.append(odeSolution.y[i+2][timeIndex])

    return x_coords, y_coords, z_coords

def prepare_all_plots():
    """Prepares plots visualising every particle at some point in time for all points in time returned by the ODE solver"""
    all_xs = []
    all_ys = []
    all_zs = []

    for i in range(odeSolution.t.size):
        x_coords, y_coords, z_coords = prepare_to_plot_instantaneous_positions(i)
        
        all_xs.append(x_coords)
        all_zs.append(z_coords)
        all_ys.append(y_coords)

    return all_xs, all_ys, all_zs

all_xs, all_ys, all_zs = prepare_all_plots()

def update_plot(time_index):
    """Updates the plot to show the galaxy at the time at time_index. time_index should be the time index, the
    returned value is the returned plot."""
    l_plot.set_data(all_xs[int(time_index)], all_ys[int(time_index)])

    return l_plot
    
axfreq = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor=axcolour)
time_slider = Slider(axfreq, 'Time index', 0.0, odeSolution.t.size -1, valinit=0, valstep=1)
time_slider.on_changed(update_plot)

def animation_frame(time_index):
    """This is run for each frame when the animation is running. It updates the plot to the relevant plot if
    the animation is running. time_index is the time index, the returned value is the returned plot."""

    if (animation_running):
        return update_plot(time_index)

anim = animation.FuncAnimation(fig, animation_frame, range(0,odeSolution.t.size -1, 1), interval = 5, repeat = True, repeat_delay = 500)
anim.event_source.stop()

def play_pause_animation(event):
    """Toggles whether or not the animation is running"""

    global animation_running
    animation_running = not animation_running

animate_button_axis = plt.axes([0.55, 0.025, 0.3, 0.04], facecolor=axcolour)
animate_button = Button(animate_button_axis, 'Play/pause animation', color=axcolour, hovercolor='0.975')
animate_button.on_clicked(play_pause_animation)

plt.show()

