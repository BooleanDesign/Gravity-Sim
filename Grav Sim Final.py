"""

Gravitation Simulation v. 0.1 Beta
Written: Nathan Diggins (BooleanDesigns)

"""
# Imports
import numpy as np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
import random as r
import os
from glib import *

"""
Getting Settings
"""
config_file = open(os.getcwd() + "\config.config", 'r+')
settings = {i.split('=')[0]: i.split('=')[1] for i in config_file.read().split('\n') if i[0] != ':'}
plt.rcParams['animation.ffmpeg_path'] = unicode(settings['writer_directory'])
G = float(settings['gravitational_constant'])
Writer = animation.writers['ffmpeg']
writer = Writer(fps=int(settings['fps']), metadata=dict(artist=settings['writer']), bitrate=int(settings['bitrate']))


def animate(i, dt=0.01):
    print i
    ax.clear()
    ax.set_xlim(cloud.axes[0], cloud.axes[1])
    ax.set_ylim(cloud.axes[2], cloud.axes[3])
    ax.set_zlim(cloud.axes[4], cloud.axes[5])
    ax.view_init(elev=60 - (i / 10.0), azim=-35 + (i / 10.0))
    cloud.update(dt=dt)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.plot(part_pos[0], part_pos[1], part_pos[2], 'o', markersize=3)
    return g


def init():
    ax.set_xlim(cloud.axes[0], cloud.axes[1])
    ax.set_ylim(cloud.axes[2], cloud.axes[3])
    ax.set_zlim(cloud.axes[4], cloud.axes[5])
    ax.view_init(elev=60, azim=-35)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.plot(part_pos[0], part_pos[1], part_pos[2], 'o', markersize=3)
    return g


cloud = Particle_Cloud([Particle(1, 1000, Vector([0, 0, 0]), Vector([0, 0, 0])),
                        Particle(1, 100, Vector([10.0, 0, 0]), Vector([0.0, 0.0, 10.0])),
                        Particle(1, 1, Vector([0, 5, 0]), Vector([10 * np.sqrt(2), 0, 0])),
                        Particle(1, 1, Vector([11, 0, 0]), Vector([np.sqrt(1000), 0, 0]))], constant_axes=True,
                       init_axes=[-15, 15, -15, 15, -15, 15])
part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
            [i.x.values[2] for i in cloud.particles]]
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax.plot(part_pos[0], part_pos[1], part_pos[2], 'o')
ani = animation.FuncAnimation(fig1, animate, frames=2000, interval=.02, init_func=init)
ani.save('test2.mp4', writer=writer)
