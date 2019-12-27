"""

Gravitation Simulation v. 0.1 Beta
Written: Nathan Diggins (BooleanDesigns)

"""
#Imports
import numpy as np
import copy
import matplotlib.pyplot as plt
from astropy import constants as const
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
import random as r
from glib import *

# Base Parameters
plt.rcParams['animation.ffmpeg_path'] = "C:\\ffmpeg\\bin\\ffmpeg.exe"
G = 6.61e-11
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15,metadata=dict(artist='Me'),bitrate=1800)


def animate(i,dt=0.01):
    print i
    ax.clear()
    ax.set_xlim(cloud.axes[0],cloud.axes[1])
    ax.set_ylim(cloud.axes[2],cloud.axes[3])
    ax.set_zlim(cloud.axes[4],cloud.axes[5])
    ax.view_init(elev=60-(i/10.0),azim=-35+(i/10.0))
    cloud.update(dt=dt)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.scatter(part_pos[0], part_pos[1], part_pos[2],s=[i.r for i in cloud.particles])
    return g

def init():
    ax.set_xlim(cloud.axes[0],cloud.axes[1])
    ax.set_ylim(cloud.axes[2],cloud.axes[3])
    ax.set_zlim(cloud.axes[4],cloud.axes[5])
    ax.view_init(elev=60,azim=-35)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.scatter(part_pos[0], part_pos[1], part_pos[2], s=[i.r for i in cloud.particles])
    return g


particles = [Particle(1,1,Vector([r.gauss(0,3) for i in range(3)]),Vector([0,0,0])) for g in range(300)]+[Particle(20,1000,Vector([0,0,0]),Vector([0,0,0]))]
print particles
cloud = Particle_Cloud(particles,
                       init_axes=[-5,5,-5,5,-5,5])
part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                       [i.x.values[2] for i in cloud.particles]]
fig1 = plt.figure()
ax = fig1.add_subplot(111,projection='3d')
ax.plot(part_pos[0],part_pos[1],part_pos[2],'o')
ani = animation.FuncAnimation(fig1,animate,frames=3000,interval=.02,init_func=init)

ani.save('test2.mp4',writer=writer)

