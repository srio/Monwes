from Codes.Beam import Beam
from Codes.OpticalElement import Optical_element
import matplotlib.pyplot as plt
import unittest
import numpy as np
from numpy.testing import assert_almost_equal
from Codes.Shape import *

do_plot = False

class IdealLensTest(unittest.TestCase):


    def test_ideal_lens_with_trace_ideal_lens(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  test_ideal_lens_with_trace_ideal_lens")

        beam=Beam()
        beam.set_flat_divergence(0.05,0.005)

        p=1.
        q=5.

        lens = Optical_element.initialiaze_as_ideal_lens(p,q)
        beam = lens.trace_ideal_lens(beam)

        beam.plot_xz()

        if do_plot:
            plt.show()

        assert_almost_equal(np.abs(beam.x).mean(), 0.0, 4)
        assert_almost_equal(np.abs(beam.z).mean(), 0.0, 4)


    def test_ideal_lens_with_trace_optical_element(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  test_ideal_lens_with_trace_optical_element")

        beam=Beam()
        beam.set_flat_divergence(0.05,0.005)

        p=1.
        q=5.

        lens = Optical_element.initialiaze_as_ideal_lens(p,q)
        beam = lens.trace_optical_element(beam)

        beam.plot_xz()
        if do_plot:
            plt.show()

        assert_almost_equal(np.abs(beam.x).mean(), 0.0, 4)
        assert_almost_equal(np.abs(beam.z).mean(), 0.0, 4)



    ##################### Doesn't work

    def test_ideal_lens_collimated_beam(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  test_ideal_lens_collimated_beam")

        beam=Beam()
        beam.set_circular_spot(20*1e-9)
        beam.set_divergences_collimated()
        beam.plot_xz()

        p=1.
        q=5.

        lens = Optical_element.initialiaze_as_ideal_lens(p,q,q,q)
        beam = lens.trace_optical_element(beam)

        beam.plot_xz()
        if do_plot:
            plt.show()

        assert_almost_equal(np.abs(beam.x).mean(), 0.0, 4)
        assert_almost_equal(np.abs(beam.z).mean(), 0.0, 4)


