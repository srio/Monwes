import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from monwes.Vector import Vector
#import h5py


class Beam(object):

    def __init__(self,N=10000):

        N = round(N)

        self.x = np.zeros(N)
        self.y = np.zeros(N)
        self.z = np.zeros(N)

        velocity = Vector(0.0000001, 1., 0.01)
        velocity.normalization()


        self.vx = np.ones(N) * velocity.x
        self.vy = np.ones(N) * velocity.y
        self.vz = np.ones(N) * velocity.z

        self.flag = np.zeros(N)

        self.N = N


    @classmethod
    def initialize_as_person(cls, N=10000):
        beam1 = Beam(round(N *0.25))
        beam1.set_point(0. + 100, 0., 20. + 100)
        beam1.set_circular_spot(5.)

        beam2 = Beam(round(N *0.25))
        beam2.set_point(0. + 100, 0., 0. + 100)
        beam2.set_rectangular_spot(20., -20., 15., 10.)

        beam = beam1.merge(beam2)

        beam3 = Beam(round(N *0.5))
        beam3.set_point(0. + 100, 0., 0. + 100)
        beam3.set_rectangular_spot(5., -5., 10., -40.)

        beam = beam.merge(beam3)

        [beam.x, beam.z] = [beam.x-np.mean(beam.x), beam.z-np.mean(beam.z)]
        [beam.x, beam.z] = [beam.x*0.05, beam.z*0.05]


        return beam


    @classmethod
    def import_from_file(cls,filename='filename'):

        filename = 'dati/Beam/' + filename + '.h5'
        f = h5py.File(filename, 'r')

        n = np.ones(1)
        f["/N"].read_direct(n)

        beam = Beam(int(n[0]))

        f["/x"].read_direct(beam.x)
        f["/y"].read_direct(beam.y)
        f["/z"].read_direct(beam.z)
        f["/vx"].read_direct(beam.vx)
        f["/vy"].read_direct(beam.vy)
        f["/vz"].read_direct(beam.vz)
        f["/flag"].read_direct(beam.flag)

        f.close()

        return beam

    def set_point(self,x,y,z):

        self.x = x + self.x
        self.y = y + self.y
        self.z = z + self.z


    def initialize_from_arrays(self, x,y,z,vx,vy,vz,flag):

        self.x = x
        self.y = y
        self.z = z

        self.vx = vx
        self.vy = vy
        self.vz = vz

        self.flag = flag

        self.N = x.size


    def duplicate(self):
        b = Beam(N=self.N)
        b.initialize_from_arrays(self.x.copy(),
                                 self.y.copy(),
                                 self.z.copy(),
                                 self.vx.copy(),
                                 self.vy.copy(),
                                 self.vz.copy(),
                                 self.flag.copy(),
                                 )
        return b

    def good_beam(self):
        indices = np.where(self.flag>=0)
        N=np.size(indices)
        beam = Beam(N)
        beam.x= self.x[indices].copy()
        beam.y= self.y[indices].copy()
        beam.z= self.z[indices].copy()
        beam.vx= self.vx[indices].copy()
        beam.vy= self.vy[indices].copy()
        beam.vz= self.vz[indices].copy()
        beam.flag= self.flag[indices].copy()

        print("Good beam indices")
        print(indices)

        return beam

    def part_of_beam(self,indices):

        N=np.size(indices)
        beam = Beam(N)
        beam.x= self.x[indices].copy()
        beam.y= self.y[indices].copy()
        beam.z= self.z[indices].copy()
        beam.vx= self.vx[indices].copy()
        beam.vy= self.vy[indices].copy()
        beam.vz= self.vz[indices].copy()
        beam.flag= self.flag[indices].copy()

        return beam




    def number_of_good_rays(self):
        return np.size(np.where(self.flag >= 0))


    def set_rectangular_spot(self,xmax,xmin,zmax,zmin):
        self.x = (np.random.random(self.N)-0.5)*(xmax-xmin)+(xmax+xmin)/2 + self.x
        self.z = (np.random.random(self.N)-0.5)*(zmax-zmin)+(zmax+zmin)/2 + self.z


    def set_circular_spot(self,r1):
        theta = (np.random.random(self.N))*2*np.pi
        r = np.sqrt(np.random.random(self.N))*r1
        self.x = r*np.cos(theta) + self.x
        self.z = r*np.sin(theta) + self.z

    def set_gaussian_spot(self,dx,dz):
        N=self.N
        self.x = dx * (np.random.randn(N))
        self.z = dz * (np.random.randn(N))


    def set_gaussian_divergence(self,dx,dz):                                                                      # gaussian velocity distribution
        N=self.N
        self.vx = dx * (np.random.randn(N))
        self.vz = dz * (np.random.randn(N))
        self.vy = np.random.random(N)
        self.vy = np.sqrt(1 - self.vx**2 - self.vz**2)


    def set_flat_divergence(self,dx,dz):                                                                         # uniform velocity distribution
        N=self.N
        self.vx = dx * (np.random.random(N) - 0.5)*2
        self.vz = dz * (np.random.random(N) - 0.5)*2
        self.vy = np.random.random(N)
        self.vy = np.sqrt(1 - self.vx**2 - self.vz**2)


    def set_flat_divergence_with_different_optical_axis(self,dx,dz):

        N = self.N

        x = self.vx
        z = self.vz
        y = self.vy


        self.vx = dx * (np.random.random(N) - 0.5)*2
        self.vz = dz * (np.random.random(N) - 0.5)*2
        self.vx = self.vx + x
        self.vz = self.vz + z


        velocity = Vector(self.vx, self.vy, self.vz)
        velocity.normalization()

        self.vx = velocity.x
        self.vy = velocity.y
        self.vz = velocity.z




    def set_divergences_collimated(self):
            self.vx = self.x * 0.0
            self.vz = self.z * 0.0
            self.vy = self.y * 0.0 + 1.


    def merge(self,beam):

        beam_out=Beam(self.N+beam.N)

        beam_out.x[0:self.N] = self.x
        beam_out.y[0:self.N] = self.y
        beam_out.z[0:self.N] = self.z
        beam_out.vx[0:self.N] = self.vx
        beam_out.vy[0:self.N] = self.vy
        beam_out.vz[0:self.N] = self.vz
        beam_out.flag[0:self.N] = self.flag

        beam_out.x[self.N:self.N+beam.N] = beam.x
        beam_out.y[self.N:self.N+beam.N] = beam.y
        beam_out.z[self.N:self.N+beam.N] = beam.z
        beam_out.vx[self.N:self.N+beam.N] = beam.vx
        beam_out.vy[self.N:self.N+beam.N] = beam.vy
        beam_out.vz[self.N:self.N+beam.N] = beam.vz
        beam_out.flag[self.N:self.N+beam.N] = beam.flag


        return beam_out

    def retrace(self,distance,resetY=False):


        t = distance / self.vy


        self.x = self.x + t * self.vx
        self.z = self.z + t * self.vz

        if resetY:
            self.y *= 0.0
        else:
            self.y = self.y + t * self.vy



    def save_on_file(self,filename):

        filename = 'dati/Beam/' + filename +'.h5'
        f = h5py.File(filename, 'w')

        f["N"] = self.N
        f["x"] = self.x
        f["y"] = self.y
        f["z"] = self.z
        f["vx"] = self.vx
        f["vy"] = self.vy
        f["vz"] = self.vz
        f["flag"] = self.flag

        f.close()

    #
    #
    #  graphics
    #
    def plot_xz(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(self.x, self.z,color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('x axis')
        plt.ylabel('z axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_xy(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(self.x, self.y,color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('x axis')
        plt.ylabel('y axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_zy(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(self.z, self.y,color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('z axis')
        plt.ylabel('y axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_yx(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(self.y, self.x,color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('y axis')
        plt.ylabel('x axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_xpzp(self,equal_axis=1, color='b', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(1e6*self.vx, 1e6*self.vz,color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('xp axis [urad]')
        plt.ylabel('zp axis [urad]')
        if equal_axis==1:
            plt.axis('equal')

    def plot_ypzp(self,equal_axis=1, color='b', marker='.', markersize=0.2):
        plt.figure()
        plt.plot(1e6*self.vy, 1e6*self.vz, color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('yp axis [urad]')
        plt.ylabel('zp axis [urad]')
        if equal_axis==1:
            plt.axis('equal')


    def plot_good_xz(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        plt.figure()
        indices = np.where(self.flag >= 0)
        plt.plot(self.x[indices], self.z[indices], color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('x axis')
        plt.ylabel('z axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_good_xy(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        indices = np.where(self.flag >= 0)
        plt.figure()
        plt.plot(self.x[indices], self.y[indices],color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('x axis')
        plt.ylabel('y axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_good_zy(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        indices = np.where(self.flag >= 0)
        plt.figure()
        plt.plot(self.z[indices], self.y[indices], color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('z axis')
        plt.ylabel('y axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_good_yx(self,equal_axis=1, color='r', marker='.', markersize=0.2):
        indices = np.where(self.flag >= 0)
        plt.figure()
        plt.plot(self.y[indices], self.x[indices], color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('y axis')
        plt.ylabel('x axis')
        if equal_axis==1:
            plt.axis('equal')

    def plot_good_xpzp(self,equal_axis=1, color='b', marker='.', markersize=0.2):
        indices = np.where(self.flag >= 0)
        plt.figure()
        plt.plot(1e6*self.vx[indices], 1e6*self.vz[indices], color=color, marker=marker, markersize=markersize, linestyle='None')
        plt.xlabel('xp axis [urad]')
        plt.ylabel('zp axis [urad]')
        if equal_axis==1:
            plt.axis('equal')


    def histogram(self):
        plt.figure()
        plt.hist(self.x,100)
        plt.title('x position for gaussian distribution')

        plt.figure()
        plt.hist(self.z,100)
        plt.title('z position for uniform distribution')


