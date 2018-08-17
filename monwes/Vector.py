import numpy as np
import matplotlib.pyplot as plt
from numpy import pi

class Vector(object):

    def __init__(self, x, y, z):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

    def duplicate(self):
        return Vector(self.x,self.y,self.z)

    def size(self):
        return self.x.size

    def rotation_x(self,alpha):
        x=self.x.copy()
        y=self.y.copy()
        z=self.z.copy()

        self.x = x
        self.y = y*np.cos(alpha)-z*np.sin(alpha)
        self.z = y*np.sin(alpha)+z*np.cos(alpha)

    def rotation_y(self,beta):
        x=self.x.copy()
        y=self.y.copy()
        z=self.z.copy()
        self.x = x*np.cos(beta)+z*np.sin(beta)
        self.y = y
        self.z = -x*np.sin(beta)+z*np.cos(beta)


    def rotation_z(self,gamma):
        x=self.x.copy()
        y=self.y.copy()
        z=self.z.copy()
        self.x = x*np.cos(gamma)-y*np.sin(gamma)
        self.y = x*np.sin(gamma)+y*np.cos(gamma)
        self.z = z



    def rotation(self,angle,axis="x"):
        """
        rotate a vector an angle alpha
        :param alpha: rotation angle in degrees (counterclockwise)
        :param axis: "x", "y" or "z"
        :return:
        """

        if axis == "x":
            self.rotation_x(angle)
        elif axis=="y":
            self.rotation_y(angle)
        elif axis=="z":
            self.rotation_z(angle)


    def surface_conic_normal(self,ccc):

        x=2*ccc[1-1]*self.x+ccc[4-1]*self.y+ccc[6-1]*self.z+ccc[7-1]
        y=2*ccc[2-1]*self.y+ccc[4-1]*self.x+ccc[5-1]*self.z+ccc[8-1]
        z=2*ccc[3-1]*self.z+ccc[5-1]*self.y+ccc[6-1]*self.x+ccc[9-1]

        return Vector(x,y,z)


    def modulus(self):
        return np.sqrt(self.x**2+self.y**2+self.z**2)


    def normalization(self):
        mod = self.modulus()
        if self.x.size != 1:
            indices = np.where(np.abs(mod)>1e-13)
            self.x[indices] = self.x[indices] / mod[indices]
            self.y[indices] = self.y[indices] / mod[indices]
            self.z[indices] = self.z[indices] / mod[indices]

        else:
            if mod > 1e-13:
                self.x = self.x / mod
                self.y = self.y / mod
                self.z = self.z / mod
            else:
                self.x = 0.
                self.y = 0.
                self.z = 0.


    def dot(self,v2):
        return np.array(self.x*v2.x+self.y*v2.y+self.z*v2.z)

    def vector_product(self,w):
        return Vector(self.y * w.z - self.z * w.y,
                      self.z * w.x - self.x * w.z,
                      self.x * w.y - self.y * w.x)

    def perpendicular_component(self,normal):
        a=-self.dot(normal)
        return Vector(
            normal.x*a,
            normal.y*a,
            normal.z*a)

    def sum(self,v2):
        return Vector(  self.x+v2.x,
                        self.y+v2.y,
                        self.z+v2.z)


    def rodrigues_formula(self,axis1,theta):

        mod = self.modulus()

        axis = axis1.duplicate()
        axis.normalization()
        vrot=Vector(self.x,self.y,self.z)
        vrot.x=self.x*np.cos(theta)+( axis.y*self.z-axis.z*self.y)*np.sin(theta)+(1-np.cos(theta))*axis.x**2*self.x
        vrot.y=self.y*np.cos(theta)+(-axis.x*self.z+axis.z*self.x)*np.sin(theta)+(1-np.cos(theta))*axis.y**2*self.y
        vrot.z=self.z*np.cos(theta)+( axis.x*self.y-axis.y*self.x)*np.sin(theta)+(1-np.cos(theta))*axis.z**2*self.z

        indices = np.where(np.abs(mod)>1e-13)
        vrot.normalization()

        if self.x.size != 1:
            vrot.x[indices] *= mod[indices]
            vrot.y[indices] *= mod[indices]
            vrot.z[indices] *= mod[indices]

        return vrot



    def info(self):
        if self.size() == 1:
            return "x: %f, y: %f, z: %f\n"%(self.x,self.y,self.z)
        else:
            txt = ""
            for i in range(self.size()):
                txt += "x: %f, y: %f, z: %f\n"%(self.x[i],self.y[i],self.z[i])
            return txt

