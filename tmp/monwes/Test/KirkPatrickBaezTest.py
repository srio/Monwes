from monwes.Beam import Beam
import unittest
from monwes.Shape import BoundaryRectangle
import numpy as np
import matplotlib.pyplot as plt
from monwes.CompoundOpticalElement import CompoundOpticalElement
from numpy.testing import assert_almost_equal

do_plot = False

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

    oe0.FDISTR = 1
    oe0.FSOUR = 0
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.005
    oe0.HDIV2 = 0.005
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 10000
    oe0.PH1 = 1000.0
    oe0.VDIV1 = 0.05
    oe0.VDIV2 = 0.05


    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam







class KirkPatrickBaezTest(unittest.TestCase):

    def test_kirk_patrick_baez(self):

        shadow_beam = shadow_source()
        beam = Beam()
        beam.initialize_from_arrays(
            shadow_beam.getshonecol(1),
            shadow_beam.getshonecol(2),
            shadow_beam.getshonecol(3),
            shadow_beam.getshonecol(4),
            shadow_beam.getshonecol(5),
            shadow_beam.getshonecol(6),
            shadow_beam.getshonecol(10),
        )



        bound1 = BoundaryRectangle(xmax=2.5, xmin=-2.5, ymax=2.5, ymin=-2.5)
        bound2 = BoundaryRectangle(xmax=1., xmin=-1., ymax=1., ymin=-1.)

        kirk_patrick_baez = CompoundOpticalElement.initialize_as_kirkpatrick_baez(p=10., q=5., separation=4., theta=89*np.pi/180, bound1=bound1, bound2=bound2)

        beam = kirk_patrick_baez.trace_compound(beam)


        beam.plot_good_xz(0)

        indices = np.where(beam.flag>0)
        assert_almost_equal(beam.x[indices], 0., 4)
        assert_almost_equal(beam.z[indices], 0., 4)

        beam.retrace(50.)

        beam.plot_good_xz()

        print(kirk_patrick_baez.info())

        if do_plot:
            plt.show()
