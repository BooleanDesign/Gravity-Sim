# Imports
import numpy as np
import copy
import matplotlib.pyplot as plt
from astropy import constants as const
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
import random as r
import time
G = 1.0
# Class Definitions
def timer(func):
    def wrapper(*args,**kwargs):
        t1 = time.time()
        x = func(*args,**kwargs)
        t2 = time.time()
        return t2-t1,x
    return wrapper
class Vector():
    """
    Defines the Vector Class for computation
    """

    def __init__(self, data):
        """
        Initiates the vector class
        :param data: <type list> delineates the values of the vector.
        """
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

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return Vector.__mul__(self, other)

    def __sub__(self, other):
        return Vector.__add__(self, -1 * other)

    def sqrt(self):
        return Vector([np.sqrt(i) for i in self.values])

    def unit(self):
        return (1 / abs(self)) * self

    def __pow__(self, power):
        return Vector([i ** power for i in self.values])

    def __abs__(self):
        return np.sqrt(sum([i ** 2 for i in self.values]))

    def __div__(self, other):
        try:
            return self.__mul__(1.0 / other)
        except:
            raise TypeError

    def cross(self, other):
        return Vector([np.linalg.det([[self.values[1], self.values[2]],
                                      [other.values[1], other.values[2]]]),
                       -1 * np.linalg.det([[self.values[0], self.values[2]],
                                           [other.values[0], other.values[2]]]),
                       np.linalg.det([[self.values[1], self.values[2]],
                                      [other.values[1], other.values[2]]])])


class particle():
    """

    Initiates the particle class, which acts as a shaped point particle with typical physical properties

    """

    def __init__(self, radius, mass, position, velocity):
        self.r = radius
        self.m = mass
        self.x = position
        self.v = velocity
        self.a = 0

    def update(self, dt, force):
        """
        This updates the position of the particle.
        :param dt: This is the time difference to calculate over
        :param force: This is the force of the particle
        :return: None
        """
        self.a = force / self.m
        self.v += dt * self.a
        self.x += dt * self.v

    def __add__(self, other):
        if type(self) == type(other):
            print (self.m * self.v).values, (other.m * other.v).values, self.m + other.m
            return particle(self.r + other.r, self.m + other.m, self.x,
                            ((self.m * self.v) + (other.m * other.v)) / (self.m + other.m))
        else:
            return self
    def __radd__(self,other):
        return self.__add__(other)
    def distance(self, other):
        return abs(self.x - other.x)

    def momentum(self):
        return self.m * self.v


class Particle_Cloud():
    """
    Defines the aggregate particle cloud. This object consists of a group of particles.
    """

    def __init__(self, particles, axes_parameter=1.1, size_parameter=0.008, constant_axes=True, init_axes=None):
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
            self.axes = [(2.0 - axes_parameter) * min([i.x.values[0] for i in self.particles]),
                         axes_parameter * max([i.x.values[0] for i in self.particles]),
                         (2.0 - axes_parameter) * min([i.x.values[1] for i in self.particles]),
                         axes_parameter * max([i.x.values[1] for i in self.particles]),
                         (2.0 - axes_parameter) * min([i.x.values[2] for i in self.particles]),
                         axes_parameter * max([i.x.values[2] for i in self.particles])]
        if (self.axes[1] - self.axes[0]) * (self.axes[3] - self.axes[2]) * (self.axes[5] - self.axes[4]) < 2.0:
            self.axes = [-2, 2, -2, 2, -2, 2]
        self.lin_mom = sum([i.m * i.v for i in self.particles])  # Defines the linear momentum vector
        self.com = (1 / sum([i.m for i in particles])) * (
            sum([i.m * i.x for i in particles]))  # Defines the center of mass
        self.ang_mom = sum([self.com.cross(i.m * i.v) for i in self.particles])
    @timer
    def update2(self, dt=0.01):
        part_sims = []
        """
        Removing merged particles
        """
        new_particles = []
        for part in self.particles:
            if part in self.particles:
                non_self_particles = [i for i in self.particles if i != part]
                part_sims = []
                for i, j in enumerate([part.distance(k) < self.size_parameter for k in non_self_particles]):
                    if j == True and non_self_particles[i] in self.particles:
                        part_sims.append(i)
                        self.particles.remove(part)
                        new_particles.append(part + non_self_particles[i])

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
        for part in parts:
            force = gravitational_force(part, parts)
            self.particles[parts.index(part)].update(dt, force)
        if self.constant_axes == True:
            try:
                self.axes = self.init_axes
            except:
                self.axes = [-2, 2, -2, 2, -2, 2]
        else:
            self.axes = [(2.0 - self.axes_param) * min([i.x.values[0] for i in self.particles]),
                         self.axes_param * max([i.x.values[0] for i in self.particles]),
                         (2.0 - self.axes_param) * min([i.x.values[1] for i in self.particles]),
                         self.axes_param * max([i.x.values[1] for i in self.particles]),
                         (2.0 - self.axes_param) * min([i.x.values[2] for i in self.particles]),
                         self.axes_param * max([i.x.values[2] for i in self.particles])]
        if (self.axes[1] - self.axes[0]) * (self.axes[3] - self.axes[2]) * (self.axes[5] - self.axes[4]) < 2.0:
            self.axes = [-2, 2, -2, 2, -2, 2]
        self.lin_mom = sum([i.m * i.v for i in self.particles])  # Defines the linear momentum vector
        self.com = (1 / sum([i.m for i in self.particles])) * (
            sum([i.m * i.x for i in self.particles]))  # Defines the center of mass
        self.ang_mom = sum([self.com.cross(i.m * i.v) for i in self.particles])
    @timer
    def update(self, dt=0.01):
        part_sims = []
        """
        Removing merged particles
        """
        new_particles = []
        removable_particles = []
        for p in range(len(self.particles)): #This allows iteration through the whole list of particles
            temp_removable_particles = []
            distances = [self.particles[p].distance(i) <= (self.particles[p].r+i.r) for i in self.particles[p+1:]]
            print distances
            if True in distances:
                temp_removable_particles = [self.particles[p]]
                temp_removable_particles += [self.particles[(p+1)+i] for i in range(len(self.particles[p+1:])) if distances[i] == True] #Checks to see if a given particle is in the distances list
            new_particles.append(sum(temp_removable_particles))
            removable_particles += temp_removable_particles
        for particle in removable_particles:
            self.particles.remove(particle)
        for particle in new_particles:
            if type(particle) != int: #Checks that the sumation did not yield a zero value.
                self.particles.append(particle)
        """
        Updating particles
        """
        parts = copy.deepcopy(self.particles)
        for part in parts:
            force = gravitational_force(part, parts)
            self.particles[parts.index(part)].update(dt, force)
        if self.constant_axes:
            try:
                self.axes = self.init_axes
            except:
                self.axes = [-2, 2, -2, 2, -2, 2]
        else:
            self.axes = [(2.0 - self.axes_param) * min([i.x.values[0] for i in self.particles]),
                         self.axes_param * max([i.x.values[0] for i in self.particles]),
                         (2.0 - self.axes_param) * min([i.x.values[1] for i in self.particles]),
                         self.axes_param * max([i.x.values[1] for i in self.particles]),
                         (2.0 - self.axes_param) * min([i.x.values[2] for i in self.particles]),
                         self.axes_param * max([i.x.values[2] for i in self.particles])]
        if (self.axes[1] - self.axes[0]) * (self.axes[3] - self.axes[2]) * (self.axes[5] - self.axes[4]) < 2.0:
            self.axes = [-2, 2, -2, 2, -2, 2]
        self.lin_mom = sum([i.m * i.v for i in self.particles])  # Defines the linear momentum vector
        self.com = (1 / sum([i.m for i in self.particles])) * (
            sum([i.m * i.x for i in self.particles]))  # Defines the center of mass
        self.ang_mom = sum([self.com.cross(i.m * i.v) for i in self.particles])


def gravitational_force(i, n):
    """
    Returns the gravitational force on the particle.
    :param i: <type particle> this is the base particle
    :param n: <type list(particle)> this is the particle cloud values
    :return: returns the force on the particle <Vector>
    """
    w = [b for b in n if i != b] #takes all the objects in b except for those in n
    f = 1.0 * i.m * float(G) * sum([j.m * ((j.x - i.x).unit()) / (abs(j.x - i.x) ** 2) for j in w]) #Returns the net force on the particle
    return f

if __name__ == '__main__':
    h = range(1,100)
    t1 = []
    t2 = []
    for number in h:
        print number
        t = []
        g = Particle_Cloud([particle(0.1,1,Vector([r.gauss(0,2) for i in range(3)]),Vector([0,0,0])) for k in range(number)],init_axes=[-1,1,-1,1,-1,1])
        for l in range(10):
            t.append(g.update()[0])
        t1.append(sum(t)/len(t))
        t = []
        g = Particle_Cloud([particle(0.1,1,Vector([r.gauss(0,2) for i in range(3)]),Vector([0,0,0])) for k in range(number)],init_axes=[-1,1,-1,1,-1,1])
        for l in range(10):
            t.append(g.update2()[0])
        t2.append(sum(t)/len(t))
    plt.plot(h,t1)
    plt.plot(h,t2,'b')
    plt.show()
