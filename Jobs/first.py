from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
import Shadow

main = "__montel__"

def shadow_source():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 1
    oe0.FSOUR = 1
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.005
    oe0.HDIV2 = 0.005
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 0
    oe0.NPOINT = 25000
    oe0.PH1 = 9137.65
    oe0.VDIV1 = 0.00075
    oe0.VDIV2 = 0.00075
    oe0.WXSOU = 0.0003
    oe0.WZSOU = 0.00015

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam

def beam():

    shadow_beam = shadow_source()

    beam = Beam(25000)
    beam.initialize_from_arrays(
        shadow_beam.getshonecol(1),
        shadow_beam.getshonecol(2),
        shadow_beam.getshonecol(3),
        shadow_beam.getshonecol(4),
        shadow_beam.getshonecol(5),
        shadow_beam.getshonecol(6),
        shadow_beam.getshonecol(10),
    )

    beam.flag *= 0.

    #beam = Beam()
    #beam.set_divergences_collimated()
    #beam.set_point(0., 0., 0.)
    beam.set_circular_spot(1e-3)


    return beam


def plot_montel(beam, equal_axis=1):
    plt.figure()
    plt.plot(beam[0].x, beam[0].z, 'r.')
    plt.plot(beam[1].x, beam[1].z, 'b.')
    plt.plot(beam[2].x, beam[2].z, 'g.')
    plt.xlabel('x axis')
    plt.ylabel('z axis')

    if equal_axis == 1:
        plt.axis('equal')

    beam[2].plot_xz(0)


    print("No reflection = %d\nOne reflection = %d\nTwo reflection = %d" %(beam[0].N, beam[1].N, beam[2].N))

if main == "__ideal_lenses__":

    beam = beam()

    beam.plot_xz(0)
    plt.title('starting beam')

    ideal_lens_1 = Optical_element.initialiaze_as_ideal_lens(p=0.4, q=0.3, fx=0.4, fz=0.4)
    ideal_lens_2_04 = Optical_element.initialiaze_as_ideal_lens(p=0.3, q=0.4, fx=0.4, fz=0.4)
    ideal_lens_2_1471 = Optical_element.initialiaze_as_ideal_lens(p=0.3, q=1.471, fx = 1.471, fz = 1.471)


    beam = ideal_lens_1.trace_optical_element(beam)

    beam_04 = beam.duplicate()
    beam_1471 = beam.duplicate()


    beam_04 = ideal_lens_2_04.trace_optical_element(beam_04)        ########### case of 0.4
    beam_1471 = ideal_lens_2_1471.trace_optical_element(beam_1471)  ########### case of 1.471

    beam_04.plot_xz(0)
    plt.title('0.4')
    beam_1471.plot_xz(0)
    plt.title('1.471')

    plt.show()

if main == "__montel__":

    beam = beam()

    beam.x *= 1e6
    beam.z *= 1e6


    beam.plot_xpzp(0)
    beam.plot_xz(0)
    plt.title('starting beam')

    beam.x *= 1e-6
    beam.z *= 1e-6

    xmax = 0.+1e-10
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.-1e-10

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_1 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.4, q=0., theta=88.281 * np.pi / 180,
                                                                     bound1=bound,
                                                                     bound2=bound, distance_of_the_screen=0.3, fp=0.4,
                                                                     fq=10000000.)

    beam_1 = montel_1.trace_montel(beam)

    beam_1 = beam_1[2]


    #######################################################################################################################
    beam_1_04 = beam_1.duplicate()
    beam_1_14 = beam_1.duplicate()

    xmax = 0.
    xmin = -100.
    ymax = 0.4
    ymin = -0.3
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_2 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=0.4, theta=88.281 * np.pi / 180,
                                                                     bound1=bound, bound2=bound, fp=1000000000., fq=0.4)

    beam_2_04 = montel_2.trace_montel(beam_1_04)

    beam_2_04 = beam_2_04[2]

    position = Vector(beam_2_04.x, beam_2_04.y, beam_2_04.z)
    velocity = Vector(beam_2_04.vx, beam_2_04.vy, beam_2_04.vz)



    beam_2_04.x *= 1e6
    beam_2_04.z *= 1e6

    beam_2_04.plot_xz(0)


    beam_2_04.x *= 1e-6
    beam_2_04.z *= 1e-6

    plt.title('0.4')

########################################################################################################################

    xmax = 0.
    xmin = -100.
    ymax = 1.471
    ymin = -0.3
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_2 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=1.471, theta=88.281 * np.pi / 180,
                                                                     bound1=bound, bound2=bound,
                                                                     fp=1000000000., fq=1.471)

    beam_2_14 = montel_2.trace_montel(beam_1)

    beam_2_14 = beam_2_14[2]



    beam_2_14.x *= 1e6
    beam_2_14.z *= 1e6

    beam_2_14.plot_xz(0)

    beam_2_14.x *= 1e-6
    beam_2_14.z *= 1e-6

    plt.title('1.471')

    plt.show()