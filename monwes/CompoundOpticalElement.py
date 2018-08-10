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

        #oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p1,0.,theta1,alpha1,"p",2*z0-q1)
        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p1, 0., theta1, alpha1, "p", 2*z0-q1)
        #oe2 = Optical_element.initialize_my_hyperboloid(p=0.,q=q1,theta=90*np.pi/180,alpha=0,wolter=1, z0=z0, distance_of_focalization=2*z0-q1)
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
    def wolter_1_for_microscope(cls, p, q, d, q1, theta1, theta2):

        ##############  ellipse     ####################################################################################
        ae = (p+q1)/2
        be = np.sqrt(p*q1)*np.cos(theta1)
        f = np.sqrt(ae**2-be**2)
        beta = np.arccos((p**2+4*f**2-q1**2)/(4*p*f))

        ccc1 = np.array([1. / be ** 2, 1. / be ** 2, 1 / ae ** 2, 0., 0., 0., 0., 0., 0., -1])

        y = - p * np.sin(beta)
        z = f - p * np.cos(beta)

        oe1 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc1)
        oe1.set_parameters(p=p, q=q1, theta=theta1)
        ##############   hyperbola  ####################################################################################
        p1 = q1 - d
        ah = (p1 - q)/2
        bh = np.sqrt(p1*q)*np.cos(theta2)
        z0 = np.sqrt(ae**2-be**2) - np.sqrt(ah**2+bh**2)

        print("z0 = %f" %z0)

        ccc2 = np.array([-1. / ah ** 2, -1. / ah ** 2, 1 / bh ** 2, 0., 0., 0., 0., 0., 2 * z0 / bh ** 2, z0 ** 2 / bh ** 2 - 1])

        oe2 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc2)
        oe2.set_parameters(p=0., q=q, theta=90*np.pi/180)

        return CompoundOpticalElement(oe_list=[oe1, oe2], oe_name="Wolter 1 for microscope")



    @classmethod
    def initialize_as_kirkpatrick_baez(cls, p, q, separation, theta, bound1, bound2):


        p1 = p - 0.5 * separation
        q1 = p - p1
        q2 = q - 0.5 * separation
        p2 = q - q2
        f1p = p1
        f1q = p+ q - p1
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
    def initialize_as_montel_parabolic(cls, p, q, theta, bound1, bound2, distance_of_the_screen=None, angle_of_mismatch=0.):

        beta = (90. - angle_of_mismatch)*np.pi/180    #### angle beetween the two mirror, if angle_of_mismatch is >0 the two mirror are closer


        oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=p, q=q, theta=theta, alpha=0., infinity_location='p', focal=q, cylindrical=1)
        oe1.set_bound(bound1)

        oe2 = oe1.duplicate()
        oe2.rotation_surface_conic(beta, 'y')
        oe2.set_bound(bound2)

        if distance_of_the_screen == None:
            distance_of_the_screen = q
        ccc = np.array([0., 0., 0., 0., 0., 0., 0., 1., 0., -distance_of_the_screen])
        screen = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
        screen.set_parameters(p, distance_of_the_screen, 0., 0., "Surface conical mirror")



        return CompoundOpticalElement(oe_list=[oe1, oe2, screen], oe_name="Montel parabolic")

    @classmethod
    def initialize_as_montel_ellipsoid(cls, p, q, theta, bound1, bound2, distance_of_the_screen=None, angle_of_mismatch=0., fp=None, fq=None):

        beta = (90.- angle_of_mismatch)*np.pi/180    #### angle beetween the two mirror


        oe1 = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p=p, q=q, theta=theta, alpha=0., cylindrical=1, fp=fp, fq=fq)
        oe1.set_bound(bound1)

        oe2 = oe1.duplicate()
        oe2.rotation_surface_conic(beta, 'y')
        oe2.set_bound(bound2)


        if distance_of_the_screen == None:
            distance_of_the_screen = q

        print(distance_of_the_screen)

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

                self.oe[i].set_parameters(theta=180.*np.pi/180, q=q1*0.)





    def compound_specification_after_screen(self, beam, i):
        if self.type == "Wolter 1":
            if i == 1:
                self.oe[i].output_frame_wolter(beam)


        if self.type == "Wolter 2":
            if i == 1:
                self.oe[i].output_frame_wolter(beam)

        if self.type == "Wolter 3":
            if self.oe[i].theta < 1e-10:
                self.oe[i].output_frame_wolter(beam)
                #todo control well this part
                x = beam.x
                z = beam.z
                vx = beam.vx
                vz = beam.vz

                beam.x = z
                beam.z = x
                beam.vx = vz
                beam.vz = vx

        #if self.type == "Wolter 1 for microscope":
        #    if i == 1:
        #        self.oe[i].output_frame_wolter(beam)



    def system_initialization(self, beam):

        if self.type == "Wolter 1 for microscope":
            print("Wolter for japanese")
            self.velocity_wolter_microscope(beam)



    def trace_compound(self,beam1):

        beam=beam1.duplicate()

        self.system_initialization(beam)

        for i in range (self.oe_number()):

            print("Iteration number %d"  %(i+1))

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

    def rotation_traslation_montel(self, beam):

        theta = self.oe[0].theta
        p = self.oe[0].p



        theta = np.pi / 2 - theta

        vector = Vector(0., 1., 0.)
        vector.rotation(-theta, 'x')


        ny = -vector.z / np.sqrt(vector.y ** 2 + vector.z ** 2)
        nz = vector.y / np.sqrt(vector.y ** 2 + vector.z ** 2)

        n = Vector(0, ny, nz)

        vrot = vector.rodrigues_formula(n, -theta)
        vrot.normalization()


        #########################################################################################################################

        position = Vector(beam.x, beam.y, beam.z)
        mod_position = position.modulus()
        velocity = Vector(beam.vx, beam.vy, beam.vz)

        position.rotation(-theta, 'x')
        velocity.rotation(-theta, 'x')

        position = position.rodrigues_formula(n, -theta)
        velocity = velocity.rodrigues_formula(n, -theta)
        velocity.normalization()

        #position.normalization()
        position.x = position.x #* mod_position
        position.y = position.y #* mod_position
        position.z = position.z #* mod_position

        [beam.x, beam.y, beam.z] = [position.x, position.y, position.z]
        [beam.vx, beam.vy, beam.vz] = [velocity.x, velocity.y, velocity.z]

        ####### translation  ###################################################################################################

        vector_point = Vector(0, p, 0)

        vector_point.rotation(-theta,  "x")
        vector_point = vector_point.rodrigues_formula(n, -theta)
        vector_point.normalization()

        beam.x = beam.x - vector_point.x * p
        beam.y = beam.y - vector_point.y * p
        beam.z = beam.z - vector_point.z * p

        return beam


    def time_comparison(self, beam1, elements):

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

    def trace_montel(self, beam):

        beam = self.rotation_traslation_montel(beam)

        # beam.plot_xz()

        origin = self.time_comparison(beam, elements=[1, 2, 3])

        indices = np.where(origin == 1)
        beam1 = beam.part_of_beam(indices)
        indices = np.where(origin == 2)
        beam2 = beam.part_of_beam(indices)
        indices = np.where(origin == 3)
        beam3 = beam.part_of_beam(indices)

        # beam1.plot_xz(0)
        # beam2.plot_xz(0)
        # beam3.plot_xz(0)

        origin0 = origin


        if beam3.N != 0:
            [beam3, t] = self.oe[2].intersection_with_optical_element(beam3)

        beam1_list = [beam1.duplicate(), Beam(), Beam()]
        beam2_list = [beam2.duplicate(), Beam(), Beam()]
        beam3_list = [beam3.duplicate(), Beam(), Beam()]


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


            beam3_list[i + 1] = beam003.merge(beam03)
            if beam3_list[i + 1].N != 0:
                [beam3_list[i + 1], t] = self.oe[2].intersection_with_optical_element(beam3_list[i + 1])


        print("Super origin")

        print(origin0)
        print(origin1)
        print(origin2)

        self.print_file_montel(beam1_list, beam2_list, beam3_list, origin0, origin1, origin2)

        #self.oe[0].output_frame_wolter(beam3_list[2])

        return beam3_list


    def print_file_montel(self,beam_1, beam_2, beam_3, origin0, origin1, origin2):

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

        indeces_x2_none = indices

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

        f = open("dati.txt", "w")

        f.write(
            "Index     ons     dos     tros         x1          y1          z1        x2        y2         z2         x3         y3          z3\n")

        for i in range(0, len(x3)):

            # print(x2[i])

            if x1[i] == np.inf:
                a1 = "     "
            else:
                a1 = ""

            if x2[i] == np.inf:
                a2 = "    "
            else:
                a2 = ""

            if x3[i] == np.inf:
                if x2[i] != np.inf:
                    a3 = "   "
                else:
                    a3 = ""
            else:
                a3 = ""

            if i < 10:
                f.write((
                                    "%d           %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e2:
                f.write((
                                    "%d          %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e3:
                f.write((
                                    "%d         %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e4:
                f.write((
                                    "%d        %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e5:
                f.write((
                                    "%d       %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e6:
                f.write((
                                    "%d      %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))
            elif i < 1e7:
                f.write((
                                    "%d     %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" % (
                            i, unos[i], dos[i], tros[i], a1, x1[i], a1, a1, y1[i], a1, a1, z1[i], a1, a2, x2[i], a2, a2,
                            y2[i], a2, a2, z2[i], a2, a3, x3[i], a3, a3, y3[i], a3, a3, z3[i], a3)))

        f.close()

        import h5py

        f = h5py.File("dati.h5",'w')

        f.attrs['NX_class'] = 'NXentry'
        f.attrs['default'] = 'montel_footprint'


        f1 = f.create_group("montel_footprint")


        f1.attrs['NX_class'] = 'NXdata'
        f1.attrs['signal'] = b'z1'
        f1.attrs['axes'] = b'x1'

        f1["i"] = range(0, len(x3))
        f1["unos"] = unos
        f1["x1"] = x1
        f1["y1"] = y1
        f1["z1"] = z1
        f1["z3"] = z3
        f.close()

        f = h5py.File("dati.h5",'r')
        X1 = f["montel_footprint/x1"].value
        f.close()

        print("Hello World")
        print(X1-x1)

    def velocity_wolter_microscope(self, beam):


        ccc = self.oe[0].ccc_object.get_coefficients()

        ae = ccc[2]**-0.5
        be = ccc[1]**-0.5


        b2 = beam.y
        b3 = beam.z
        beam.y = b3
        beam.z = b2
        beam.z = - beam.z + np.sqrt(ae ** 2 - be ** 2)



        p = self.oe[0].p
        q = self.oe[0].q
        print("ccc")
        print(ccc,ccc[2]**-0.5,ccc[0]**-0.5)
        ae = ccc[2]**-0.5
        be = 1/np.sqrt(ccc[0])
        f = np.sqrt(ae**2-be**2)
        print("p=%f, q=%f, ae=%f, be=%f f=%f" %(p,q,ae,be,f))

        beta = np.arccos((p**2+4*f**2-q**2)/(4*p*f))
        y = p * np.sin(beta)
        z = f - p * np.cos(beta)

        v = Vector(0., y/p, (z-f)/p)
        v.normalization()
        v0 = Vector(0., 0., -1.)
        v0.normalization()
        alpha = np.arccos(v.dot(v0))

        velocity = Vector(beam.vx, beam.vz, -beam.vy)
        velocity.rotation(-alpha, 'x')

        beam.vx = velocity.x
        beam.vy = velocity.y
        beam.vz = velocity.z



        #position = Vector(beam.x, beam.z, -beam.y)
        #position.rotation(-alpha, 'x')

        #beam.x = position.x
        #beam.y = position.y
        #beam.z = position.z

        t = (-v.x * beam.x - v.y * beam.y - v.z * beam.z) / (
                v.x * beam.vx + v.y * beam.vy + v.z * beam.vz)
        beam.x += beam.vx * t
        beam.y += beam.vy * t
        beam.z += beam.vz * t


        self.oe[0].set_parameters(p=0., q=0., theta=90*np.pi/180)



    def trace_wolter_japanese(self,beam):

        self.velocity_wolter_microscope(beam)
        self.oe[0].intersection_with_optical_element(beam)
        self.oe[0].output_direction_from_optical_element(beam)
        self.oe[1].intersection_with_optical_element(beam)
        self.oe[1].output_direction_from_optical_element(beam)

        ccc = self.oe[1].ccc_object.get_coefficients()

        ah = (-ccc[0])**-0.5
        bh = ccc[2]**-0.5
        z0 = ccc[8]*bh**2/2

        t = (-z0+np.sqrt(ah**2+bh**2)-beam.z) / beam.vz

        print("ah = %f, bh = %f" %(ah, bh))

        beam.x += beam.vx * t
        beam.y += beam.vy * t
        beam.z += beam.vz * t

        self.oe[0].output_frame_wolter(beam)

        return beam