from Codes.Beam import Beam
from Codes.Shape import BoundaryRectangle
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Codes.CompoundOpticalElement import CompoundOpticalElement

do_plot = True
main = "__main__"


if main == "__main__":


    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_wolter2")

    p = 100.  ##### if p=100 the trace_good_ray goes crazy
    beam1 = Beam.initialize_as_person(10000)
    # beam1 = Beam(100000)
    # beam1.set_circular_spot(1.)
    # beam1.set_rectangular_spot(5 / 2 * 1e-5, -5 / 2 * 1e-5, 5 / 2 * 1e-5, -5 / 2 * 1e-5)
    beam1.x *= 100.
    beam1.z *= 100.
    beam1.set_point(p, 0., p)

    op_ax = Beam(1)
    op_ax.set_point(p, 0., p)

    beam = op_ax.merge(beam1)
    beam.set_divergences_collimated()
    beam.plot_xz(0)

    p = 20000.
    q = 30.
    z0 = 5.
    focal = 2 * z0 + q

    wolter2 = CompoundOpticalElement.initialiaze_as_wolter_2(p1=p, q1=q, z0=z0)

    beam = wolter2.trace_good_rays(beam)

    beam.plot_good_xz()

    beam.retrace(10.)
    beam.plot_good_xz()
    plt.title("test_wolter2_good_rays")

    print(beam.flag)

    if do_plot:
        plt.show()
