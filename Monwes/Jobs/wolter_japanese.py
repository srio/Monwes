from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt

main = '__wolter_japanese_v0.2__'

def shadow_source():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
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

    oe0.FDISTR = 3
    oe0.FSOUR = 1
    oe0.F_PHOT = 0
    oe0.HDIV1 = 5.9999998e-05
    oe0.HDIV2 = 5.9999998e-05
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 25000
    oe0.PH1 = 1000.0
    oe0.SIGDIX = 4.6699999e-05
    oe0.SIGDIZ = 0.00025499999
    oe0.VDIV1 = 0.00028499999
    oe0.VDIV2 = 0.00028499999
    oe0.WXSOU = 0.04
    oe0.WZSOU = 0.002

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
    beam.x *= 1e-2
    beam.z *= 1e-2

    return beam

if main == '__wolter_japanese_v0.1__':

    beam = beam()


    beam.x *= 1e6
    beam.z *= 1e6

    #beam.plot_xz(0)
    #beam.plot_xpzp(0)

    beam.x *= 1e-6
    beam.z *= 1e-6


    wolter_japanese = CompoundOpticalElement.wolter_1_for_microscope(p=13.4, q=0.3, q1=0.67041707, d=0.1082882, theta1=88.8*np.pi/180, theta2=89.*np.pi/180)


    beam = wolter_japanese.trace_wolter_japanese(beam)

    wolter_japanese.oe[0].output_frame_wolter(beam)

    beam.x *= 1e6
    beam.z *= 1e6

    beam.plot_xz(0)
    beam.plot_xpzp(0)

    print(beam.vx, beam.vy, beam.vz)

    plt.show()


if main == '__wolter_japanese_v0.2__':

    beam = beam()
    wolter_japanese = CompoundOpticalElement.wolter_1_for_microscope(p=13.4, q=0.3, q1=0.67041707, d=0.1082882, theta1=88.8*np.pi/180, theta2=89.*np.pi/180)

########## Initialization  of the beam  ################################################################################

    ccc = wolter_japanese.oe[0].ccc_object.get_coefficients()

    ae = ccc[2] ** -0.5
    be = ccc[1] ** -0.5
    f = np.sqrt(ae ** 2 - be ** 2)

    b2 = beam.y
    b3 = beam.z
    beam.y = b3
    beam.z = b2
    beam.z = - beam.z + f

    p = wolter_japanese.oe[0].p
    q = wolter_japanese.oe[0].q

    beta = np.arccos((p ** 2 + 4 * f ** 2 - q ** 2) / (4 * p * f))
    y = - p * np.sin(beta)
    z = f - p * np.cos(beta)

    v = Vector(0., y, z - f)
    v.normalization()
    v0 = Vector(0., 0., -1.)
    v0.normalization()
    alpha = np.arccos(v.dot(v0))


    t = (-v.x * beam.x - v.y * beam.y - v.z * beam.z) / (
            v.x * beam.vx + v.y * beam.vy + v.z * beam.vz)
    beam.x += beam.vx * t
    beam.y += beam.vy * t
    beam.z += beam.vz * t


    velocity = Vector(beam.vx, beam.vz, -beam.vy)
    velocity.rotation(-alpha, 'x')

    beam.vx = velocity.x
    beam.vy = velocity.y
    beam.vz = velocity.z

    wolter_japanese.oe[0].set_parameters(p=0., q=0., theta=90 * np.pi / 180)

########################################################################################################################

    wolter_japanese.oe[0].intersection_with_optical_element(beam)
    wolter_japanese.oe[0].output_direction_from_optical_element(beam)
    wolter_japanese.oe[1].intersection_with_optical_element(beam)
    wolter_japanese.oe[1].output_direction_from_optical_element(beam)


    ccc = wolter_japanese.oe[1].ccc_object.get_coefficients()

    ah = (-ccc[0]) ** -0.5
    bh = ccc[2] ** -0.5
    z0 = ccc[8] * bh ** 2 / 2

    #t = (-z0 + np.sqrt(ah ** 2 + bh ** 2) - beam.z) / beam.vz
    t = - beam.y / beam.vy

    print("ah = %f, bh = %f" % (ah, bh))

    beam.x += beam.vx * t
    beam.y += beam.vy * t
    beam.z += beam.vz * t

    beam.plot_xz()

    print(wolter_japanese.info())


    plt.show()
