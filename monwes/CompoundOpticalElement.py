from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
import numpy as np
import matplotlib.pyplot as plt
from monwes.Vector import Vector
import time


class CompoundOpticalElement(object):

    def __init__(self,oe_list=[],oe_name=""):
        self.oe = oe_list
        self.type = oe_name

    def append_oe(self,oe):
        self.oe.append(oe)

    def oe_number(self):
        return len(self.oe)

    def reset_oe_list(self):
        self.oe = []

    def set_type(self,name):
        self.type = name


    @classmethod
    def initialiaze_as_wolter_1(cls,p1,q1,z0):
        theta1 = 0.
        alpha1 = 0.
        print(q1)
        print(2*z0)
        print("dof=%f" %(2*z0-q1))

        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p1, 0., theta1, alpha1, "p", 2*z0-q1)
        oe2 = Optical_element.initialize_my_hyperboloid(p=0., q=q1, theta=90 * np.pi / 180, alpha=0, wolter=1, z0=z0,distance_of_focalization=2*z0-q1)

        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Wolter 1")


    @classmethod
    def initialiaze_as_wolter_1_with_two_parameters(cls,p1, R, theta):

        cp1 = -2 * R / np.tan(theta)
        cp2 = 2 * R * np.tan(theta)
        cp = max(cp1, cp2)
        f = cp / 4
        print("focal=%f" % (f))

        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=p1, q=f, theta=0., alpha=0.,infinity_location="p")

        s1 = R / np.tan(2 * theta)
        s2 = R / np.tan(4 * theta)
        c = (s1 - s2) / 2
        z0 = f + c


        b1 = np.sqrt(
            0.5 * c ** 2 + 0.5 * R ** 2 + 0.5 * R ** 4 / cp ** 2 - R ** 2 * z0 / cp + 0.5 * z0 ** 2 - 0.5 / cp ** 2 * np.sqrt(
                (
                            -c ** 2 * cp ** 2 - cp ** 2 * R ** 2 - R ** 4 + 2 * cp * R ** 2 * z0 - cp ** 2 * z0 ** 2) ** 2 - 4 * cp ** 2 * (
                            c ** 2 * R ** 4 - 2 * c ** 2 * cp * R ** 2 * z0 + c ** 2 * cp ** 2 * z0 ** 2)))
        b2 = np.sqrt(
            0.5 * c ** 2 + 0.5 * R ** 2 + 0.5 * R ** 4 / cp ** 2 - R ** 2 * z0 / cp + 0.5 * z0 ** 2 + 0.5 / cp ** 2 * np.sqrt(
                (
                            -c ** 2 * cp ** 2 - cp ** 2 * R ** 2 - R ** 4 + 2 * cp * R ** 2 * z0 - cp ** 2 * z0 ** 2) ** 2 - 4 * cp ** 2 * (
                            c ** 2 * R ** 4 - 2 * c ** 2 * cp * R ** 2 * z0 + c ** 2 * cp ** 2 * z0 ** 2)))
        b = min(b1, b2)
        a = np.sqrt(c ** 2 - b ** 2)

        ccc = np.array(
            [-1 / a ** 2, -1 / a ** 2, 1 / b ** 2, 0., 0., 0., 0., 0., -2 * z0 / b ** 2, z0 ** 2 / b ** 2 - 1])
        oe2 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
        oe2.set_parameters(p=0., q=z0+c, theta=90*np.pi/180, alpha=0., type="My hyperbolic mirror")
        #oe2.type = "My hyperbolic mirror"
        #oe2.p = 0.
        #oe2.q = z0 + c
        #oe2.theta = 90 * np.pi / 180
        #oe2.alpha = 0.


        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Wolter 1")


    @classmethod
    def initialiaze_as_wolter_2(cls,p1,q1,z0):
        #q1 = - q1
        focal = q1+2*z0
        print("focal=%f" %(focal))
        theta1 = 0.
        alpha1 = 0.

        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p1,0.,theta1,alpha1,"p", focal)
        oe2 = Optical_element.initialize_my_hyperboloid(p=0. ,q=-(focal-2*z0), theta=90*np.pi/180, alpha=0, wolter=2, z0=z0, distance_of_focalization=focal)


        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Wolter 2")



    @classmethod
    def initialiaze_as_wolter_12(cls,p1,q1,focal_parabola,Rmin):


        focal = focal_parabola
        d = q1 - focal_parabola
        z0 = focal_parabola + d/2
        print("focal=%f" %(focal))
        theta1 = 0.
        alpha1 = 0.


        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p1,0.,theta1,alpha1,"p", focal)
        ccc = oe1.ccc_object.get_coefficients()
        cp = -ccc[8]
        print("R=%f, d=%f, cp=%f, z0=%f" %(Rmin,d,cp,z0))
        #b1 = np.sqrt(0.125*d**2+Rmin+2*Rmin**4/cp**2-2*Rmin**2*z0/cp+0.5*z0**2-0.125/cp**2*np.sqrt((-cp**2*d**2-8*cp**2*Rmin-16*Rmin**4+16*cp*Rmin**2*z0-4*cp*z0**2)**2-16*cp**2*(4*d**2*Rmin**4-4*cp*d**2*Rmin**2*z0+cp**2*d**2*z0**2)))
        #b2 = np.sqrt(0.125*d**2+Rmin+2*Rmin**4/cp**2-2*Rmin**2*z0/cp+0.5*z0**2+0.125/cp**2*np.sqrt((-cp**2*d**2-8*cp**2*Rmin-16*Rmin**4+16*cp*Rmin**2*z0-4*cp*z0**2)**2-16*cp**2*(4*d**2*Rmin**4-4*cp*d**2*Rmin**2*z0+cp**2*d**2*z0**2)))
        p1 = -cp ** 2 * d ** 2 - 8 * cp ** 2 * Rmin - 16 * Rmin ** 4 + 16 * cp * Rmin ** 2 * z0 - 4 * cp ** 2 * z0 ** 2
        p1 = p1**2
        p2 = 16 * cp ** 2 * (4 * d ** 2 * Rmin ** 4 - 4 * cp * d ** 2 * Rmin ** 2 * z0 + cp ** 2 * d ** 2 * z0 ** 2)
        sp = 0.125/cp**2*np.sqrt(p1-p2)
        sp0 = 0.125*d**2+Rmin+2*Rmin**4/cp**2-2*Rmin**2*z0/cp+0.5*z0**2
        b = np.sqrt(sp0-sp)
        a = np.sqrt(d**2/4-b**2)

        print("a=%f, b=%f" %(a,b))


        #oe2 = Optical_element.initialize_my_hyperboloid(p=0. ,q=-(focal-2*z0), theta=90*np.pi/180, alpha=0, wolter=1.1, z0=z0, distance_of_focalization=focal)

        cc = np.array([-1/a**2, -1/a**2, 1/b**2, 0., 0., 0., 0., 0., -2*z0/b**2, (z0/b)**2-1])
        oe2 = Optical_element.initialize_as_surface_conic_from_coefficients(cc)
        oe2.type = "My hyperbolic mirror"
        oe2.set_parameters(p=0., q=q1, theta=90.*np.pi/180, alpha=0.)

        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Wolter 1.2")


    @classmethod
    def initialize_as_wolter_3(cls, p, q, distance_between_the_foci):
        f=-q

        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=p, q=f, theta=0., alpha=0., infinity_location="p")

        #c = z0+np.abs(f)
        c = distance_between_the_foci/2
        z0 = np.abs(c)-np.abs(f)
        b = c+100
        a = np.sqrt((b**2-c**2))
        ccc = np.array([1/a**2, 1/a**2, 1/b**2, 0., 0., 0., 0., 0., -2*z0/b**2, z0**2/b**2-1])

        oe2 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
        oe2.set_parameters(p=0., q=z0+z0+np.abs(q), theta=90*np.pi/180)


        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Wolter 3")



    @classmethod
    def initialize_as_kirkpatrick_baez(cls, p, q, separation, theta, bound1=None, bound2=None):


        p1 = p - 0.5 * separation
        q1 = p - p1
        q2 = q - 0.5 * separation
        p2 = q - q2
        f1p = p1
        f1q = p + q - p1
        f2q = q2
        f2p = p + q - q2


        oe1 = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p= f1p, q= f1q, theta= theta, alpha=0., cylindrical=1)
        #oe1.bound = bound1
        oe1.set_bound(bound1)
        oe1.p = p1
        oe1.q = q1

        oe2 = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p= f2p, q= f2q, theta= theta, alpha=90.*np.pi/180, cylindrical=1)
        #oe2.bound = bound2
        oe2.set_bound(bound2)
        oe2.p = p2
        oe2.q = q2

        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Kirkpatrick Baez")



    @classmethod
    def initialize_as_kirkpatrick_baez_parabolic(cls, p, q, separation, theta, bound1=None, bound2=None, infinity_location='p'):


        p1 = p - 0.5 * separation
        q1 = p - p1
        q2 = q - 0.5 * separation
        p2 = q - q2
        f1p = p1
        f1q = p + q - p1
        f2q = q2
        f2p = p + q - q2


        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p= f1p, q= f1q, theta= theta, alpha=0., cylindrical=1, infinity_location=infinity_location)
        #oe1.bound = bound1
        oe1.set_bound(bound1)
        oe1.p = p1
        oe1.q = q1

        oe2 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p= f2p, q= f2q, theta= theta, alpha=90.*np.pi/180, cylindrical=1, infinity_location=infinity_location)
        #oe2.bound = bound2
        oe2.set_bound(bound2)
        oe2.p = p2
        oe2.q = q2

        print(oe1.ccc_object.get_coefficients())
        print(oe2.ccc_object.get_coefficients())

        return CompoundOpticalElement(oe_list=[oe1,oe2],oe_name="Kirkpatrick Baez")




    @classmethod
    def initialize_as_montel_parabolic(cls, p, q, theta_z, theta_x=None, bound1=None, bound2=None, distance_of_the_screen=None, angle_of_mismatch=0., infinity_location='q'):

        beta = np.pi / 2 + angle_of_mismatch            #### angle beetween the two mirror, if angle_of_mismatch is <0 the two mirror are closer
        if theta_x == None:
            theta_x = theta_z

        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=p, q=q, theta=theta_z, alpha=0.*np.pi/180, infinity_location=infinity_location, cylindrical=1.)
        oe2 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=p, q=q, theta=theta_x, alpha=0.*np.pi/180, infinity_location=infinity_location, cylindrical=1.)
        oe1.set_bound(bound1)

        oe2.rotation_surface_conic(beta, 'y')
        oe2.set_bound(bound2)



        if distance_of_the_screen == None:
            distance_of_the_screen = q
        ccc = np.array([0., 0., 0., 0., 0., 0., 0., 1., 0., -distance_of_the_screen])
        screen = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
        screen.set_parameters(p, distance_of_the_screen, 0., 0., "Surface conical mirror")

        return CompoundOpticalElement(oe_list=[oe1, oe2, screen], oe_name="Montel parabolic")

    @classmethod
    def initialize_as_montel_ellipsoid(cls, p, q, theta_z, theta_x=None, bound1=None, bound2=None, angle_of_mismatch=0., fp=None, fq=None):

        beta = np.pi / 2 + angle_of_mismatch            #### angle beetween the two mirror, if angle_of_mismatch is <0 the two mirror are closer
        if theta_x == None:
            theta_x = theta_z

        oe1 = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p=p, q=q, theta=theta_z, alpha=0., cylindrical=1, fp=fp, fq=fq)
        oe2 = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p=p, q=q, theta=theta_x, alpha=0., cylindrical=1, fp=fp, fq=fq)
        oe1.set_bound(bound1)


        oe2.rotation_surface_conic(beta, 'y')
        oe2.set_bound(bound2)

        distance_of_the_screen = q


        ccc = np.array([0., 0., 0., 0., 0., 0., 0., 1., 0., -distance_of_the_screen])
        screen = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
        screen.set_parameters(p, q, 0., 0., "Surface conical mirror")


        return CompoundOpticalElement(oe_list=[oe1, oe2, screen], oe_name="Montel ellipsoid")

    def compound_specification_after_oe(self, i):
        if self.type == "Wolter 1":
            if i == 0:
                self.oe[i].set_parameters(p=None, q=0., theta=90.*np.pi/180)
            elif i == 1:
                self.oe[i].set_parameters(p=None, q=None, theta=0.)


        if self.type == "Wolter 1.2":
            if i == 0:
                self.oe[i].set_parameters(p=None, q=0., theta=90.*np.pi/180)
            elif i == 1:
                self.oe[i].set_parameters(p=None, q=None, theta=0.)



        if self.type == "Wolter 2":
            if i == 0:
                self.oe[i].set_parameters(p=None, q=0., theta=90.*np.pi/180)
            elif i == 1:
                self.oe[i].set_parameters(p=None, q=None, theta=0.)

        if self.type == "Wolter 3":
            if np.abs(self.oe[i].theta) < 1e-10:
                self.oe[i].set_parameters(p=None, q=0., theta=90.*np.pi/180)
            else:
                self.oe[i].set_parameters(p=None, q=None, theta=0.)


        if self.type == "Wolter 1 for microscope":
            if i == 1:
                print("oh yeah")
                ccc = self.oe[i].ccc_object.get_coefficients()
                ah = (-ccc[0]) ** -0.5
                bh = ccc[2] ** -0.5
                z0 = ccc[8] * bh ** 2 / 2

                q1 = z0 - np.sqrt(ah ** 2 + bh ** 2)

                print("The coordinate of the hyperbola focus is %g" %q1)

                self.oe[i].set_parameters(theta=90.*np.pi/180, q=q1*0.)



    def compound_specification_after_screen(self, beam, i):
        if self.type == "Wolter 1":
            if i == 1:

                vx = np.mean(beam.vx)
                vy = np.mean(beam.vy)
                vz = np.mean(beam.vz)

                v = Vector(vx, vy, vz)

                self.output_frame(beam, v, mode=1.)
                print(np.mean(beam.vx), np.mean(beam.vy), np.mean(beam.vz))




        if self.type == "Wolter 2":
            if i == 1:

                vx = np.mean(beam.vx)
                vy = np.mean(beam.vy)
                vz = np.mean(beam.vz)

                v = Vector(vx, vy, vz)

                self.output_frame(beam, v, mode=1.)


        if self.type == "Wolter 3":
            if self.oe[i].theta < 1e-10:

                vx = np.mean(beam.vx)
                vy = np.mean(beam.vy)
                vz = np.mean(beam.vz)

                v = Vector(vx, vy, vz)

                self.output_frame(beam, v, mode=1.)


    def system_initialization(self, beam):

        if self.type == "Wolter 1 for microscope":
            print("Wolter for japanese")
            self.velocity_wolter_microscope(beam)



    def trace_compound(self,beam1):

        beam=beam1.duplicate()

        self.system_initialization(beam)

        for i in range (self.oe_number()):

            print("\n\nIteration number %d"  %(i+1))

            if self.oe[i] is not None:

                self.oe[i].effect_of_optical_element(beam)
                self.compound_specification_after_oe(i=i)
                self.oe[i].effect_of_the_screen(beam)
                self.compound_specification_after_screen(beam = beam, i=i)

        return beam


    def info(self):

        txt = ("\nThe optical element of the %s system are:\n" %(self.type))

        for i in range (self.oe_number()):
            txt += ("\nThe %d optical element:\n\n" %(i+1))
            txt += self.oe[i].info()
        return txt



    def trace_good_rays(self, beam1):

        beam11=beam1.duplicate()
        beam = beam1.duplicate()


        self.oe[0].rotation_to_the_optical_element(beam11)
        self.oe[0].translation_to_the_optical_element(beam11)

        b1=beam11.duplicate()
        b2=beam11.duplicate()
        [b1, t1] = self.oe[0].intersection_with_optical_element(b1)
        [b2, t2] = self.oe[1].intersection_with_optical_element(b2)


        indices = np.where(beam.flag>=0)
        beam.flag[indices] = beam.flag[indices] + 1

        if self.type == "Wolter 1":
            indices = np.where (t1>=t2)
        elif self.type == "Wolter 2":
            indices = np.where (t1<=t2)

        beam.flag[indices] = -1*beam.flag[indices]
        print(beam.flag)


        print("Trace indices")
        indices = np.where(beam.flag>=0)
        print(indices)
        #beam.plot_good_xz(0)

        beam = beam.good_beam()

        beam.plot_good_xz()
        plt.title("Good initial rays")


        l = beam.number_of_good_rays()
        print(l)

        if l >0:
            beam = self.trace_compound(beam)
        else:
            print(">>>>>>NO GOOD RAYS")


        print("Number of good rays=%f" %(beam.number_of_good_rays()))

        return beam



    def time_comparison(self, beam1, elements):

        #
        # take as input a beam and optical elements and compute the travel time that every ray
        # has to do to reach the optical elements.
        # The output is a vector that indicate which is the closest optical element for every ray
        #

        origin = np.ones(beam1.N)
        tf = 1e35 * np.ones(beam1.N)
        tau = np.ones(3)*6

        for i in range (0, len(elements)):


            beam = beam1.duplicate()

            [beam, t] = self.oe[elements[i]-1].intersection_with_optical_element(beam)
            #t = abs(t)
            tau[i] = np.mean(t)

            indices = np.where(beam.flag<0)


            t[indices] = 1e30


            tf = np.minimum(t, tf)
            indices = np.where(t == tf)
            origin[indices] = elements[i]

        return origin



    def get_optical_axis_in(self, mode=0):

        theta_grazing_z = np.pi/ 2 - self.oe[0].theta
        theta_grazing_x = np.pi/ 2 - self.oe[1].theta

        if mode == 0:
            vz = np.tan(theta_grazing_z) / np.sqrt(1 + np.tan(theta_grazing_z) ** 2)
            vx = np.tan(theta_grazing_x) / np.sqrt(1 + np.tan(theta_grazing_x) ** 2)


        if mode == 1:
            vz = np.tan(theta_grazing_z) / np.sqrt(1 + np.tan(theta_grazing_z) ** 2 + np.tan(theta_grazing_x) ** 2)
            vx = np.tan(theta_grazing_x) / np.sqrt(1 + np.tan(theta_grazing_z) ** 2 + np.tan(theta_grazing_x) ** 2)


        return Vector(vx, np.sqrt(1 - vx**2 - vz**2), -vz)



    def get_optical_axis_angles(self, v, mode=0):

        v = v.duplicate()

        alpha = np.arctan(v.x/v.y)
        v.rotation(alpha, 'z')
        beta = np.arctan(-v.z/v.y)
        v.rotation(beta, 'x')

        return alpha, beta

    def rotation(self, beam, angle, axis):
        velocity = Vector(beam.vx, beam.vy, beam.vz)
        position = Vector(beam.x, beam.y, beam.z)

        velocity.rotation(angle, axis)
        position.rotation(angle, axis)

        beam.x = position.x
        beam.y = position.y
        beam.z = position.z


        beam.vx = velocity.x
        beam.vy = velocity.y
        beam.vz = velocity.z


    def translation_before_the_system(self, beam, v, hitting_point):

        beam.x -= self.oe[0].p * v.x - hitting_point.x
        beam.y -= self.oe[0].p * v.y - hitting_point.y
        beam.z -= self.oe[0].p * v.z - hitting_point.z



    def apply_specular_reflections(self, beam, name_file, print_footprint):

        # This is the core of the tracing algoritm


        # This first part compute the travelling time before the two mirror and the image plane
        # and so divide the beam in three:
        # beam1: are the rays that wil first hit oe1
        # beam2: are the rays that will first hit oe2
        # beam3: are the rays that doesn't hit any mirror

        origin = self.time_comparison(beam, elements=[1, 2, 3])
        indices = np.where(origin == 1)
        beam1 = beam.part_of_beam(indices)
        indices = np.where(origin == 2)
        beam2 = beam.part_of_beam(indices)
        indices = np.where(origin == 3)
        beam3 = beam.part_of_beam(indices)

        origin0 = origin.copy()


        if beam3.N != 0:
            [beam3, t] = self.oe[2].intersection_with_optical_element(beam3)

        beam1_list = [beam1.duplicate(), Beam(), Beam()]
        beam2_list = [beam2.duplicate(), Beam(), Beam()]
        beam3_list = [beam3.duplicate(), Beam(), Beam()]




        # This second part compute calculate the travel time of the ray, after the first reflection,
        # to reach the others elements, and put the rays in those beam
        # beam1_list: are the rays that hit oe1
        # beam2_list: are the rays that hit oe2
        # beam3_list: are the rays that reach the image plane


        origin1 = [1,2]
        origin2 = [1,2]

        for i in range(0, 2):

            if beam1_list[i].N != 0:
                [beam1_list[i], t] = self.oe[0].intersection_with_optical_element(beam1_list[i])
                self.oe[0].output_direction_from_optical_element(beam1_list[i])

                origin = self.time_comparison(beam1_list[i], [2, 3])
                origin1[i] = origin
                indices = np.where(origin == 2)
                beam2_list[i + 1] = beam1_list[i].part_of_beam(indices)
                indices = np.where(origin == 3)
                beam03 = beam1_list[i].part_of_beam(indices)
            else:
                beam2_list[i+1] = Beam(0)
                beam03 = Beam(0)

            if beam2_list[i].N != 0:
                [beam2_list[i], t] = self.oe[1].intersection_with_optical_element(beam2_list[i])
                self.oe[1].output_direction_from_optical_element(beam2_list[i])

                origin = self.time_comparison(beam2_list[i], [1, 3])
                origin2[i] = origin
                indices = np.where(origin == 1)
                beam1_list[i + 1] = beam2_list[i].part_of_beam(indices)
                indices = np.where(origin == 3)
                beam003 = beam2_list[i].part_of_beam(indices)
            else:
                beam1_list[i+1] = Beam(0)
                beam003 = Beam(0)


            beam3_list[i + 1] = beam03.merge(beam003)



        print("Resuming")
        print(beam1_list[0].N, beam1_list[1].N, beam1_list[2].N)
        print(beam2_list[0].N, beam2_list[1].N, beam2_list[2].N)
        print(beam3_list[0].N, beam3_list[1].N, beam3_list[2].N)

        self.print_file_montel(beam_1=beam1_list, beam_2=beam2_list, beam_3=beam3_list, origin0=origin0,
                               origin1=origin1, origin2=origin2, name_file=name_file, print_footprint=print_footprint)

        return  beam1_list, beam2_list, beam3_list


    def get_optical_axis_out(self, mode):

        if mode == 0:
            theta_grazing = np.pi/ 2 - self.oe[0].theta
            vx = np.tan(theta_grazing) / np.sqrt(1 + np.tan(theta_grazing) ** 2)
            return Vector(-vx, np.sqrt(1 - 2 * vx **2), vx)

        if mode == 1:
            theta_grazing = np.pi/ 2 - self.oe[0].theta
            vx = np.tan(theta_grazing) / np.sqrt(1 + 2 * np.tan(theta_grazing) ** 2)
            return Vector(-vx, np.sqrt(1 - 2 * vx **2), vx)


    def translation_after_the_system(self, beam):

        beam.y -= self.oe[1].q

        t = - beam.y / beam.vy
        beam.x += beam.vx * t
        beam.y += beam.vy * t
        beam.z += beam.vz * t


    def output_frame(self, beam, v_out, mode):

        #
        # change the frame with one that depend from the the vector v_out
        #

        alpha, beta = self.get_optical_axis_angles(v_out, mode)
        self.rotation(beam, alpha, 'z')
        self.rotation(beam, beta,  'x')
        self.translation_after_the_system(beam)


    def input_frame(self, beam, v_in, mode, hitting_point):
        #
        # This method put the center of the frame with the center of montel
        #

        alpha, beta = self.get_optical_axis_angles(v_in, mode)

        self.rotation(beam, -beta,  'x')
        self.rotation(beam, -alpha, 'z')
        self.translation_before_the_system(beam, v_in, hitting_point)



    def trace_montel(self, beam, name_file=None, mode=0, p=None, q=None, theta_z=None, theta_x=None, hitting_point=Vector(0., 0., 0.), output_frame=0., print_footprint=1):


        #
        # p: is the source distance from the montel, q is the image plane distance from the montel, theta_z correspond to the incidence angle with the first yz-mirror
        # theta_x with the xy-mirror
        # hitting point are the coordinate where the beam hit the montel system
        # output fram = 0: image_plane that correspond to the two reflection beam
        # output fram = 1: image_plane that correspond to the no reflection beam
        #
        self.oe[0].set_parameters(p=p, q=q, theta=theta_z)
        self.oe[1].set_parameters(p=p, q=q, theta=theta_x)

        v_in = self.get_optical_axis_in(mode)

        self.input_frame(beam, v_in, mode, hitting_point)


        beam1, beam2, beam3 = self.apply_specular_reflections(beam, name_file, print_footprint)

        if output_frame == 0:
            v_out = self.get_optical_axis_out(mode)
        elif output_frame == 1:
            v_out = v_in

        self.output_frame(beam3[0], v_out, mode)
        self.output_frame(beam3[1], v_out, mode)
        self.output_frame(beam3[2], v_out, mode)



        return beam1, beam2, beam3




    def print_file_montel(self,beam_1, beam_2, beam_3, origin0, origin1, origin2, name_file, print_footprint):

        unos = origin0

        dos = origin0 * 0.
        indices = np.where(unos == 3)
        dos[indices] = dos[indices] * 0.
        indices = np.where(unos == 1)
        dos[indices] = origin1[0]
        indices = np.where(unos == 2)
        dos[indices] = origin2[0]

        tros = origin0 * 0.
        indices = np.where(dos == 3)
        tros[indices] = tros[indices] * 0.
        indices = np.where(dos == 1)
        tros[indices] = origin1[1]
        indices = np.where(dos == 2)
        tros[indices] = origin2[1]

        x1 = np.ones(len(unos)) * 0.
        y1 = np.ones(len(unos)) * 0.
        z1 = np.ones(len(unos)) * 0.

        indices = np.where(unos == 1)
        x1[indices] = beam_1[0].x
        y1[indices] = beam_1[0].y
        z1[indices] = beam_1[0].z
        indices = np.where(unos == 2)
        x1[indices] = beam_2[0].x
        y1[indices] = beam_2[0].y
        z1[indices] = beam_2[0].z
        indices = np.where(unos == 3)
        x1[indices] = beam_3[0].x
        y1[indices] = beam_3[0].y
        z1[indices] = beam_3[0].z
        indices = np.where(unos == 0)
        x1[indices] = np.inf
        y1[indices] = np.inf
        z1[indices] = np.inf

        x2 = np.ones(len(unos)) * 0.
        y2 = np.ones(len(unos)) * 0.
        z2 = np.ones(len(unos)) * 0.

        indices = np.where(dos == 1)
        x2[indices] = beam_1[1].x
        y2[indices] = beam_1[1].y
        z2[indices] = beam_1[1].z
        indices = np.where(dos == 2)
        x2[indices] = beam_2[1].x
        y2[indices] = beam_2[1].y
        z2[indices] = beam_2[1].z
        indices = np.where(dos == 3)
        x2[indices] = beam_3[1].x
        y2[indices] = beam_3[1].y
        z2[indices] = beam_3[1].z
        indices = np.where(dos == 0)
        x2[indices] = np.inf
        y2[indices] = np.inf
        z2[indices] = np.inf


        x3 = np.ones(len(unos)) * 0.
        y3 = np.ones(len(unos)) * 0.
        z3 = np.ones(len(unos)) * 0.

        indices = np.where(tros == 1)
        x3[indices] = beam_1[2].x
        y3[indices] = beam_1[2].y
        z3[indices] = beam_1[2].z
        indices = np.where(tros == 2)
        x3[indices] = beam_2[2].x
        y3[indices] = beam_2[2].y
        z3[indices] = beam_2[2].z
        indices = np.where(tros == 3)
        x3[indices] = beam_3[2].x
        y3[indices] = beam_3[2].y
        z3[indices] = beam_3[2].z
        indices = np.where(tros == 0)
        x3[indices] = np.inf
        y3[indices] = np.inf
        z3[indices] = np.inf


        if name_file is not None:

            import h5py

            f = h5py.File(name_file + '.h5','w')

            #f.attrs['NX_class'] = 'NXentry'
            #f.attrs['default'] = 'montel_footprint'

            #f = h5py.File(date, 'a')


            f1 = f.create_group("montel_footprint")

            #f1.attrs['NX_class'] = 'NXdata'
            #f1.attrs['signal'] = b'z1'
            #f1.attrs['axes'] = b'x1'

            f1["Number of rays"] = len(x3)
            f1["i"] = range(0, len(x3))
            f1["unos"] = unos
            f1["dos"] = dos
            f1["tros"] = tros
            f1["x1"] = x1
            f1["y1"] = y1
            f1["z1"] = z1
            f1["x2"] = x2
            f1["y2"] = y2
            f1["z2"] = z2
            f1["x3"] = x3
            f1["y3"] = y3
            f1["z3"] = z3
            #f.close()

        indices = np.where(tros==3)

        index0 = indices                      #### index of all the good rays

        indices2 = np.where(dos[indices]==1)
        index2 = index0[0][indices2]         ### this are good rays that have intersect befor with oe2 and after with oe1
        indices22 = indices2

        indices2 = np.where(dos[indices]==2)
        index1 = index0[0][indices2]         ### this are good rays that have intersect befor with oe1 and after with oe2
        indices11 = indices2

        on_do = 2*np.ones(len(index0[0]))
        on_do[indices2] -= 1

        xoe1 = np.ones(len(index0[0]))
        yoe1 = np.ones(len(index0[0]))
        zoe1 = np.ones(len(index0[0]))

        xoe1[indices11] = x1[index1]
        xoe1[indices22] = x2[index2]

        yoe1[indices11] = y1[index1]
        yoe1[indices22] = y2[index2]

        zoe1[indices11] = z1[index1]
        zoe1[indices22] = z2[index2]

        ######################################3

        xoe2 = np.ones(len(index0[0]))
        yoe2 = np.ones(len(index0[0]))
        zoe2 = np.ones(len(index0[0]))

        xoe2[indices11] = x2[index1]
        xoe2[indices22] = x1[index2]

        yoe2[indices11] = y2[index1]
        yoe2[indices22] = y1[index2]

        zoe2[indices11] = z2[index1]
        zoe2[indices22] = z1[index2]



        #f = h5py.File("dati.h5",'a')
        #f = h5py.File(date, 'a')

        if name_file is not None:
            f1 = f.create_group("montel_good_rays")

            f1["Number of rays"] = len(xoe1)
            f1["index"] = index0[0]
            f1["on_do"] = on_do
            f1["xoe1"] = xoe1
            f1["yoe1"] = yoe1
            f1["yoe2"] = yoe2
            f1["zoe2"] = zoe2

            f.close()

        if print_footprint == 1:
            self.print_footprin(on_do=on_do, xoe1=xoe1, yoe1=yoe1, zoe2=zoe2, yoe2=yoe2)



    def print_footprin(self, on_do, xoe1, yoe1, zoe2, yoe2):

        indices1 = np.where(on_do == 1)
        indices2 = np.where(on_do == 2)

        plt.figure()
        plt.plot(yoe1[indices1]*1e3, xoe1[indices1]*1e3, 'r.', markersize=0.7)
        plt.plot(yoe1[indices2]*1e3, xoe1[indices2]*1e3, 'b.', markersize=0.7)
        plt.title("footprint oe1")
        plt.xlabel("y[mm]")
        plt.ylabel("x[mm]")

        plt.figure()
        plt.plot(yoe2[indices1]*1e3, zoe2[indices1]*1e3, 'r.', markersize=0.7)
        plt.plot(yoe2[indices2]*1e3, zoe2[indices2]*1e3, 'b.', markersize=0.7)
        plt.title("footprint oe2")
        plt.xlabel("y[mm]")
        plt.ylabel("z[mm]")

