import numpy as np
import copy
import matplotlib.pyplot as plt
from astropy import constants as const
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
import random as r

plt.rcParams['animation.ffmpeg_path'] = "C:\\ffmpeg\\bin\\ffmpeg.exe"
G = 1.0
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15,metadata=dict(artist='Me'),bitrate=1800)

class Vector():
    """
    Defines the Vector Class for computation
    """

    def __init__(self, data):
        self.values = data
        self.dimensions = len(data)

    def __mul__(self, other):
        if type(other) == type(self) and self.dimensions == other.dimensions:
            return sum([self.values[i] * other.values[i] for i in range(self.dimensions)])
        else:
            try:
                return Vector([other * i for i in self.values])
            except:
                raise ValueError("Vector multiplication failed")

    def __add__(self, other):
        if type(other) != type(self):
            return self
        else:
            return Vector([self.values[i] + other.values[i] for i in range(self.dimensions)])
    def __radd__(self,other):
        return self.__add__(other)

    def __rmul__(self, other):
        return Vector.__mul__(self, other)

    def __sub__(self, other):
        return Vector.__add__(self, -1 * other)

    def sqrt(self):
        return Vector([np.sqrt(i) for i in self.values])

    def unit(self):
        return (1/abs(self))*self

    def __pow__(self, power):
        return Vector([i ** power for i in self.values])

    def __abs__(self):
        return np.sqrt(sum([i ** 2 for i in self.values]))

    def __div__(self,other):
        try:
            return self.__mul__(1.0/other)
        except:
            raise TypeError
    def cross(self,other):
        return Vector([np.linalg.det([[self.values[1],self.values[2]],
                             [other.values[1],other.values[2]]]),
                       -1*np.linalg.det([[self.values[0],self.values[2]],
                                      [other.values[0],other.values[2]]]),
                       np.linalg.det([[self.values[1],self.values[2]],
                                      [other.values[1],other.values[2]]])])

class particle():
    def __init__(self, radius, mass, position, velocity):
        self.r = radius
        self.m = mass
        self.x = position
        self.v = velocity
        self.a = 0

    def update(self,dt,force):
        """
        This updates the position of the particle.
        :param dt: This is the time difference to calculate over
        :param force: This is the force of the particle
        :return: None
        """
        self.a = force/self.m
        self.v += dt*self.a
        self.x += dt*self.v

    def __add__(self,other):
        print (self.m*self.v).values, (other.m*other.v).values, self.m+other.m
        return particle(self.r+other.r,self.m+other.m,self.x,((self.m*self.v)+(other.m*other.v))/(self.m+other.m))

    def distance(self,other):
        return abs(self.x-other.x)
    def momentum(self):
        return self.m*self.v

class particle_cloud():
    """
    Defines the aggregate particle cloud
    """
    def __init__(self,particles,axes_parameter=1.1,size_parameter=0.008,constant_axes = True,init_axes = None):
        """
        Initiates the class with list parameter of particles
        :param particles: particle list
        """
        self.particles = particles
        self.size_parameter = size_parameter
        self.axes_param = axes_parameter
        self.constant_axes = constant_axes
        self.init_axes = init_axes

        if constant_axes == True:
            try:
                self.axes = init_axes
            except:
                self.axes = [-2, 2, -2, 2, -2, 2]
        else:
            self.axes = [(2.0-axes_parameter)*min([i.x.values[0] for i in self.particles]),axes_parameter*max([i.x.values[0] for i in self.particles]),
                         (2.0 - axes_parameter)*min([i.x.values[1] for i in self.particles]),axes_parameter*max([i.x.values[1] for i in self.particles]),
                         (2.0 - axes_parameter)*min([i.x.values[2] for i in self.particles]),axes_parameter*max([i.x.values[2] for i in self.particles])]
        if (self.axes[1]-self.axes[0])*(self.axes[3]-self.axes[2])*(self.axes[5]-self.axes[4]) < 2.0:
            self.axes = [-2,2,-2,2,-2,2]
        self.lin_mom = sum([i.m*i.v for i in self.particles]) #Defines the linear momentum vector
        self.com = (1/sum([i.m for i in particles]))*(sum([i.m*i.x for i in particles])) #Defines the center of mass
        self.ang_mom = sum([self.com.cross(i.m*i.v) for i in self.particles])

    def update(self,dt = 0.01):
        part_sims = []
        """
        Removing merged particles
        """
        new_particles = []
        for particle in self.particles:
            if particle in self.particles:
                non_self_particles = [i for i in self.particles if i != particle]
                part_sims = []
                for i,j in enumerate([particle.distance(k) < self.size_parameter for k in non_self_particles]):
                    if j == True and non_self_particles[i] in self.particles:
                        part_sims.append(i)
                        self.particles.remove(particle)
                        new_particles.append(particle + non_self_particles[i])

            else:
                pass
        for removed_particle in part_sims:
            self.particles.remove(non_self_particles[removed_particle])
        for i in new_particles:
            self.particles.append(i)

        """
        Updating particles
        """
        parts = copy.deepcopy(self.particles)
        for particle in parts:
            force = gravitational_force(particle,parts)
            self.particles[parts.index(particle)].update(dt,force)
        if self.constant_axes == True:
            try:
                self.axes = self.init_axes
            except:
                self.axes = [-2, 2, -2, 2, -2, 2]
        else:
            self.axes = [(2.0-self.axes_parameter)*min([i.x.values[0] for i in self.particles]),self.axes_parameter*max([i.x.values[0] for i in self.particles]),
                         (2.0 - self.axes_parameter)*min([i.x.values[1] for i in self.particles]),self.axes_parameter*max([i.x.values[1] for i in self.particles]),
                         (2.0 - self.axes_parameter)*min([i.x.values[2] for i in self.particles]),self.axes_parameter*max([i.x.values[2] for i in self.particles])]
        if (self.axes[1]-self.axes[0])*(self.axes[3]-self.axes[2])*(self.axes[5]-self.axes[4]) < 2.0:
            self.axes = [-2,2,-2,2,-2,2]
        self.lin_mom = sum([i.m*i.v for i in self.particles]) #Defines the linear momentum vector
        self.com = (1/sum([i.m for i in self.particles]))*(sum([i.m*i.x for i in self.particles])) #Defines the center of mass
        self.ang_mom = sum([self.com.cross(i.m*i.v) for i in self.particles])

def gravitational_force(i,n):
    w = [b for b in n if i != b]
    f = 1.0*i.m*float(G)*sum([j.m*((j.x-i.x).unit())/(abs(j.x-i.x)**2) for j in w])
    return f


def animate(i,dt=0.0001):
    print i
    ax.clear()
    ax.set_xlim(cloud.axes[0],cloud.axes[1])
    ax.set_ylim(cloud.axes[2],cloud.axes[3])
    ax.set_zlim(cloud.axes[4],cloud.axes[5])
    ax.view_init(elev=60-(i/10.0),azim=-35+(i/10.0))
    cloud.update(dt=dt)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.plot(part_pos[0], part_pos[1], part_pos[2], 'o',markersize=3)
    return g

def init():
    ax.set_xlim(cloud.axes[0],cloud.axes[1])
    ax.set_ylim(cloud.axes[2],cloud.axes[3])
    ax.set_zlim(cloud.axes[4],cloud.axes[5])
    ax.view_init(elev=60,azim=-35)
    part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                [i.x.values[2] for i in cloud.particles]]
    g = ax.plot(part_pos[0], part_pos[1], part_pos[2], 'o',markersize=3)
    return g




cloud = particle_cloud([particle(1,1000,Vector([0,0,0]),Vector([0,0,0]))]+[particle(1,1,Vector([r.gauss(0,0.6) for i in range(3)]),Vector([r.gauss(1,0.3) for i in range(3)])) for k in range(100)],constant_axes=True,
                       init_axes=[-5,5,-5,5,-5,5])
part_pos = [[i.x.values[0] for i in cloud.particles], [i.x.values[1] for i in cloud.particles],
                       [i.x.values[2] for i in cloud.particles]]
fig1 = plt.figure()
ax = fig1.add_subplot(111,projection='3d')
ax.plot(part_pos[0],part_pos[1],part_pos[2],'o')
ani = animation.FuncAnimation(fig1,animate,frames=1000,interval=.02,init_func=init)
ani.save('test2.mp4',writer=writer)

