import numpy as np
import matplotlib.pyplot as plt
from monwes.Vector import Vector
from monwes.SurfaceConic import SurfaceConic



class Optical_element(object):

    def __init__(self,p=0. ,q=0. ,theta=0.0,alpha=0.0):
        self.p = p
        self.q = q
        self.theta = theta
        self.alpha = alpha
        self.type="None"
        self.bound = None

        #
        # IF ELEMENT IS native
        self.ccc_object = None # conic
        self.focal = None
        self.fx = None
        self.fz = None


    def duplicate(self):

        oe = Optical_element()

        oe.p = self.p
        oe.q = self.q
        oe.theta =self.theta
        oe.alpha =self.alpha
        oe.type = self.type
        if self.bound == None:
            oe.bound = None
        else:
            oe.bound = self.bound.duplicate()
        #
        # IF ELEMENT IS native
        oe.ccc_object = self.ccc_object.duplicate() # conic
        oe.focal = self.focal
        oe.fx = self.fx
        oe.fz = self.fz

        return oe

    def set_parameters(self,p=None,q=None,theta=None,alpha=None, type="None"):

        if p is not None:
            self.p = p
        if q is not None:
            self.q = q
        if theta is not None:
            self.theta = theta
        if alpha is not None:
            self.alpha = alpha
        if type is not "None":
            self.type = type

    def set_bound(self,bound):
        self.bound = bound


    #
    # native
    #


    @classmethod
    def initialiaze_as_ideal_lens(clsc, p, q, fx=None, fz=None):

        oe=Optical_element(p,q,0.0,0.0)
        oe.type="Ideal lens"
        if fx is None:\
            oe.fx = p*q/(p+q)
        else:
            oe.fx = fx

        if fz is None:
            oe.fz = p*q/(p+q)
        else:
            oe.fz = fz

        return oe


    #
    # conic
    #

    @classmethod
    def initialize_from_coefficients(cls):
        if np.array(cls.ccc_object).size != 10:
            raise Exception("Invalid coefficients (dimension must be 10)")
        else:
            cls.ccc_object = SurfaceConic.initialize_from_coefficients(cls)
            cls.type = "Surface conical mirror"


    @classmethod
    def initialize_as_surface_conic_plane(cls,p,q,theta=0.0,alpha=0.0):
        oe = Optical_element(p,q,theta,alpha)
        oe.type = "Surface conical mirror"
        oe.ccc_object = SurfaceConic(np.array([0,0,0,0,0,0,0,0,-1.,0]))
        return oe


    @classmethod
    def initialize_my_hyperboloid(cls,p,q,theta=0.0,alpha=0.0, wolter=None, z0=0., distance_of_focalization=0.):
        oe=Optical_element(p,q,theta,alpha)
        oe.type = "My hyperbolic mirror"
        a=q/np.sqrt(2)
        b = 1
        if wolter ==1:
            a = abs(z0 - distance_of_focalization) / np.sqrt(2)
        if wolter == 2:
            a = abs(z0-distance_of_focalization)/np.sqrt(2)
        if wolter == 1.1:
            a = distance_of_focalization/np.sqrt(2)

        print("z0=%f, distance_of_focalization=%f, a*sqrt(2)=%f" %(z0,distance_of_focalization,a*np.sqrt(2)))

        oe.ccc_object = SurfaceConic(np.array([-1, -1, 1, 0, 0, 0, 0., 0., -2 * z0, z0 ** 2 - a ** 2.]))
        return oe


    #
    # initializers from focal distances
    #

    @classmethod
    def initialize_as_surface_conic_sphere_from_focal_distances(cls, p, q, theta=0., alpha=0., cylindrical=0, cylangle=0.0,
                                                  switch_convexity=0):
        oe=Optical_element(p,q,theta,alpha)
        oe.type="Surface conical mirror"
        oe.ccc_object = SurfaceConic()
        oe.ccc_object.set_sphere_from_focal_distances(p, q, np.pi/2-theta)
        if cylindrical:
            oe.ccc_object.set_cylindrical(cylangle)
        if switch_convexity:
            oe.ccc_object.switch_convexity()
        return oe


    @classmethod
    def initialize_as_surface_conic_ellipsoid_from_focal_distances(cls, p, q, theta=0., alpha=0., cylindrical=0, cylangle=0.0,
                                                     switch_convexity=0, fp=None, fq=None):
        oe=Optical_element(p,q,theta,alpha)
        oe.type="Surface conical mirror"
        oe.ccc_object = SurfaceConic()
        if fp is None:
            fp = p
        if fq is None:
            fq=q
        oe.ccc_object.set_ellipsoid_from_focal_distances(fp, fq, np.pi/2-theta)
        if cylindrical:
            oe.ccc_object.set_cylindrical(cylangle)
        if switch_convexity:
            oe.ccc_object.switch_convexity()
        return oe


    @classmethod
    def initialize_as_surface_conic_paraboloid_from_focal_distances(cls, p, q, theta=0., alpha=0,  infinity_location="q", focal=None, cylindrical=0, cylangle=0.0,
                                                      switch_convexity=0):
        #print("p is %d" %(p))
        oe=Optical_element(p,q,theta,alpha)
        oe.type="Surface conical mirror"
        oe.focal=focal
        oe.ccc_object = SurfaceConic()
        if focal is None:
            oe.ccc_object.set_paraboloid_from_focal_distance(p, q, np.pi/2-theta, infinity_location)
        else:
            oe.ccc_object.set_paraboloid_from_focal_distance(p, focal, np.pi/2-theta, infinity_location)
        if cylindrical:
            oe.ccc_object.set_cylindrical(cylangle)
        if switch_convexity:
            oe.ccc_object.switch_convexity()

        return oe


    @classmethod
    def initialize_as_surface_conic_hyperboloid_from_focal_distances(cls, p, q, theta=0., alpha=0., cylindrical=0, cylangle=0.0,
                                                       switch_convexity=0):
        oe=Optical_element(p,q,theta,alpha)
        oe.type="Surface conical mirror"
        oe.ccc_object = SurfaceConic()
        oe.ccc_object.set_hyperboloid_from_focal_distances(p, q, np.pi/2-theta)
        if cylindrical:
            oe.ccc_object.set_cylindrical(cylangle)
        if switch_convexity:
            oe.ccc_object.switch_convexity()
        return oe

    @classmethod
    def initialize_as_surface_conic_from_coefficients(cls, ccc):
        oe = Optical_element()
        oe.type = "Surface conical mirror"
        oe.ccc_object = SurfaceConic(ccc)
        return oe



    #
    # methods to "trace"
    #


    def trace_optical_element(self, beam1):


        beam=beam1.duplicate()

        self.effect_of_optical_element(beam)

        beam.plot_yx(0)
        plt.title("footprint %s" %(self.type))

        self.effect_of_the_screen(beam)

        return beam



    def effect_of_optical_element(self,beam):

        self.rotation_to_the_optical_element(beam)
        self.translation_to_the_optical_element(beam)
        [beam, t]=self.intersection_with_optical_element(beam)
        self.output_direction_from_optical_element(beam)


    def effect_of_the_screen(self,beam):

        self.rotation_to_the_screen(beam)
        self.translation_to_the_screen(beam)
        if np.abs(self.q) > 1e-13:
            self.intersection_with_the_screen(beam)


    def intersection_with_optical_element(self, beam):
        if self.type == "Ideal lens":
            [beam, t] =self._intersection_with_plane_mirror(beam)
        elif self.type == "Surface conical mirror":
            [beam, t] =self._intersection_with_surface_conic(beam)
        elif self.type =="My hyperbolic mirror":
            [beam, t] =self._intersection_with_my_hyperbolic_mirror(beam)


        return [beam, t]


    def intersection_with_the_screen(self,beam):
        t = -beam.y / beam.vy

        beam.x = beam.x+beam.vx*t
        beam.y = beam.y+beam.vy*t
        beam.z = beam.z+beam.vz*t




    def rotation_to_the_optical_element(self, beam):

        position = Vector(beam.x,beam.y,beam.z)
        velocity = Vector(beam.vx,beam.vy,beam.vz)
        position.rotation(self.alpha,"y")
        position.rotation(-(np.pi/2-self.theta),"x")
        velocity.rotation(self.alpha,"y")
        velocity.rotation(-(np.pi/2-self.theta),"x")
        [beam.x,beam.y,beam.z] = [position.x,position.y,position.z]
        [beam.vx,beam.vy,beam.vz] = [velocity.x,velocity.y,velocity.z]


    def rotation_to_the_screen(self,beam):

        position = Vector(beam.x,beam.y,beam.z)
        velocity = Vector(beam.vx,beam.vy,beam.vz)

        if self.type == "Ideal lens":
            position.rotation((np.pi/2-self.theta),"x")
            velocity.rotation((np.pi/2-self.theta),"x")
        else:
            position.rotation(-(np.pi/2-self.theta),"x")
            velocity.rotation(-(np.pi/2-self.theta),"x")
        [beam.x,beam.y,beam.z] = [position.x,position.y,position.z]
        [beam.vx,beam.vy,beam.vz] = [velocity.x,velocity.y,velocity.z]



    def translation_to_the_optical_element(self, beam):
        vector_point=Vector(0,self.p,0)
        vector_point.rotation(self.alpha,"y")
        vector_point.rotation(-(np.pi/2-self.theta),"x")

        beam.x=beam.x-vector_point.x
        beam.y=beam.y-vector_point.y
        beam.z=beam.z-vector_point.z


    def translation_to_the_screen(self,beam):
        beam.y=beam.y-self.q


    def mirror_output_direction(self,beam):

        if self.type == "Surface conical mirror":
            normal_conic = self.ccc_object.get_normal(np.array([beam.x, beam.y, beam.z]))
        elif self.type == "My hyperbolic mirror":
            normal_conic = self.ccc_object.get_normal(np.array([beam.x, beam.y, beam.z]))
        elif self.type == "Surface conical mirror 2":
            normal_conic = self.ccc_object.get_normal(np.array([beam.x, beam.y, beam.z]))

        normal = Vector(normal_conic[0,:],normal_conic[1,:], normal_conic[2,:])
        velocity = Vector(beam.vx, beam.vy, beam.vz)
        vperp = velocity.perpendicular_component(normal)
        v2 = velocity.sum(vperp)
        v2 = v2.sum(vperp)
        [beam.vx, beam.vy, beam.vz] = [v2.x, v2.y, v2.z]




    def lens_output_direction(self,beam):

        gamma = np.arctan( beam.x/self.fx)
        alpha = np.arctan( -beam.y/self.fz)


        velocity = Vector(beam.vx, beam.vy, beam.vz)
        velocity.rotation(gamma, "y")
        velocity.rotation(alpha, "x")

        [beam.vx, beam.vy, beam.vz] = [velocity.x, velocity.y, velocity.z]




    def output_direction_from_optical_element(self, beam):

        if self.type == "Ideal lens":
            self.lens_output_direction(beam)
        else:
            self.mirror_output_direction(beam)


    def rectangular_bound(self,bound):
        self.bound=bound


    #
    # methods to intersection (all private)
    #

    def _intersection_with_plane_mirror(self, beam):

        indices = beam.flag >=0
        beam.flag[indices] = beam.flag[indices] + 1
        counter = beam.flag[indices][0]

        t=-beam.z/beam.vz
        beam.x = beam.x+beam.vx*t
        beam.y = beam.y+beam.vy*t
        beam.z = beam.z+beam.vz*t

        indices = beam.flag >=0
        beam.flag[indices] = beam.flag[indices] + 1

        if self.bound != None:
            position_x = beam.x.copy()
            indices = np.where(position_x < self.bound.xmin)
            position_x[indices] = 0
            indices = np.where(position_x > self.bound.xmax)
            position_x[indices] = 0
            indices = np.where(position_x == 0)
            beam.flag[indices] = -1 * counter

            position_y = beam.y.copy()
            indices = np.where(position_y < self.bound.ymin)
            position_y[indices] = 0
            indices = np.where(position_y > self.bound.ymax)
            position_y[indices] = 0
            indices = np.where(position_y == 0)
            beam.flag[indices] = -1 * counter


        return [beam, t]






    def _intersection_with_surface_conic(self, beam):



        indices = beam.flag >=0.
        beam.flag[indices] = beam.flag[indices] + 1
        counter = beam.flag[indices][0]



        t_source = 0.0
        [t1, t2, flag] = self.ccc_object.calculate_intercept(np.array([beam.x, beam.y, beam.z]),
                                                        np.array([beam.vx, beam.vy, beam.vz]))



        #if np.abs(np.mean(t1-t_source) >= np.abs(np.mean(t2-t_source))):
        #    print("sx")
        #    t=t1
        #else:
        #    print("sx")
        #    t=t2




        a = np.mean(t1)
        b = np.mean(t2)

        if np.mean(t1) > 0 and np.mean(t2) > 0:
            if a < b:
                t = t1
            else:
                t = t2
        elif np.mean(t1) < 0 and np.mean(t2) > 0:
           t = t2
        elif np.mean(t1) > 0 and np.mean(t2) < 0:
            t = t1
        elif np.abs(np.mean(t1)) < np.abs(np.mean(t2)):
            t = t1
        else:
            t = t2


        beam.counter = 1

        beam.x = beam.x + beam.vx * t
        beam.y = beam.y + beam.vy * t
        beam.z = beam.z + beam.vz * t

        if self.bound != None:

            position_x=beam.x.copy()
            indices = np.where(position_x < self.bound.xmin - 1e-9 * 1.)
            position_x[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_x > self.bound.xmax + 1e-9 * 1.)
            position_x[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_x == 0.0012598731458 + 1e36)
            beam.flag[indices] = -1*counter

            position_y=beam.y.copy()
            indices = np.where(position_y < self.bound.ymin - 1e-9 * 1.)
            position_y[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_y > self.bound.ymax + 1e-9 * 1.)
            position_y[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_y == 0.0012598731458 + 1e36)
            beam.flag[indices] = -1*counter


            position_z=beam.z.copy()
            indices = np.where(position_z < self.bound.zmin - 1e-9 * 1.)
            position_z[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_z > self.bound.zmax + 1e-9 * 1.)
            position_z[indices] = 0.0012598731458 + 1e36
            indices = np.where(position_z == 0.0012598731458 + 1e36)
            beam.flag[indices] = -1*counter


        return [beam, t]




    def _intersection_with_my_hyperbolic_mirror(self, beam):

        indices = beam.flag >=0
        beam.flag[indices] = beam.flag[indices] + 1
        counter = beam.flag[indices][0]


        ccc=self.ccc_object.get_coefficients()


        c0 = ccc[0]
        c1 = ccc[2]
        c2 = ccc[8]
        c3 = ccc[9]

        a = c1*beam.vz**2 + c0*(beam.vx**2 + beam.vy**2)
        b = c0*beam.x*beam.vx + c0*beam.y*beam.vy + c1*beam.z*beam.vz + c2*beam.vz/2
        c = c0*beam.x**2 + c0*beam.y**2 + c1*beam.z**2 + c2*beam.z + c3

#######################################################################################################################




        t1 = (-b - np.sqrt(b ** 2 - a * c)) / a
        t2 = (-b + np.sqrt(b ** 2 - a * c)) / a


        if  np.mean(t1)*np.mean(t2)>1:
            if np.abs(np.mean(t1)) <= np.abs(np.mean(t2)):
                t = t1
            else:
                t = t2
        elif np.mean(t1) > 0:
            t=t1
        else:
            t=t2

        beam.x = beam.x + beam.vx * t
        beam.y = beam.y + beam.vy * t
        beam.z = beam.z + beam.vz * t


        if self.bound != None:
            position_x=beam.x.copy()
            indices = np.where(position_x < self.bound.xmin)
            position_x[indices] = 0
            indices = np.where(position_x > self.bound.xmax)
            position_x[indices] = 0
            indices = np.where(position_x == 0)
            beam.flag[indices] = -1*counter

            position_y=beam.y.copy()
            indices = np.where(position_y < self.bound.ymin)
            position_y[indices] = 0
            indices = np.where(position_y > self.bound.ymax)
            position_y[indices] = 0
            indices = np.where(position_y == 0)
            beam.flag[indices] = -1*counter


        print("\n%s:   t1=%f, t2=%f" %(self.type, np.mean(t1), np.mean(t2)))

        return [beam, t]





    def trace_ideal_lens(self,beam1):

        beam=beam1.duplicate()


        t = self.p / beam.vy
        beam.x = beam.x + beam.vx * t
        beam.y = beam.y + beam.vy * t
        beam.z = beam.z + beam.vz * t

        gamma = np.arctan( beam.x/self.fx)
        alpha = np.arctan( -beam.z/self.fz)


        velocity = Vector(beam.vx, beam.vy, beam.vz)

        velocity.rotation(gamma, "z")
        velocity.rotation(alpha, "x")

        [beam.vx, beam.vy, beam.vz] = [velocity.x, velocity.y, velocity.z]


        #self.q=0
        t = self.q / beam.vy
        beam.x = beam.x + beam.vx * t
        beam.y = beam.y + beam.vy * t
        beam.z = beam.z + beam.vz * t

        beam.y *= 0.

        return beam

    def info(self):

        txt =""

        if self.type == "None":
            txt += ("No Optical element")
        if self.type == "Plane mirror":
            txt += ("%s\np=%f, q=%f, theta=%f, alpha=%f" %(self.type, self.p, self.q, self.theta, self.alpha))
        if self.type == "Spherical mirror":
            txt += ("%s\np=%f, q=%f, theta=%f, alpha=%f, R=%f" %(self.type, self.p, self.q, self.theta, self.alpha, self.R))
        if self.type == "Ideal lens":
            txt += ("%s\np=%f, q=%f, fx=%f, fz=%f" %(self.type, self.p, self.q, self.fx, self.fz))
        if self.type == "Surface conical mirror" or self.type == "My hyperbolic mirror":
            txt = ("%s\np=%f, q=%f, theta=%f, alpha=%f\n" %(self.type, self.p, self.q, self.theta, self.alpha))
            for i in range(10):
                if i != 9:
                    txt += ("c%d=%f, " %(i,self.ccc_object.get_coefficients()[i]))
                else:
                    txt += ("c%d=%f" %(i,self.ccc_object.get_coefficients()[i]))

        if self.bound is not None:
            txt += "\nThe bounds are\n"
            txt += self.bound.info()
        return txt


    def rotation_surface_conic(self, alpha, axis = 'x'):

        self.ccc_object.rotation_surface_conic(alpha, axis)


    def translation_surface_conic(self, x0, axis='x'):

        self.ccc_object.translation_surface_conic(x0, axis)