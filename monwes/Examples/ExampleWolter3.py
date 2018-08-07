from Codes.Beam import Beam
from Codes.Shape import BoundaryRectangle
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Codes.CompoundOpticalElement import CompoundOpticalElement


do_plot = True
main = "__main__"


if main == "__main__":

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_wolter3")

    p=50.
    beam1 = Beam.initialize_as_person()
    beam1.set_point(p, 0., p)
    #beam1.set_rectangular_spot(5 / 2 * 1e-5, -5 / 2 * 1e-5, 5 / 2 * 1e-5, -5 / 2 * 1e-5)

    op_ax = Beam (1)
    op_ax.set_point(p, 0., p)

    beam = op_ax.merge(beam1)
    beam.set_divergences_collimated()

    beam.plot_xz()

    distance_between_the_foci = 10.

    wolter3 = CompoundOpticalElement.initialize_as_wolter_3(20., 5., distance_between_the_foci)

    print(wolter3.oe[0].ccc_object.get_coefficients())
    print(wolter3.oe[1].ccc_object.get_coefficients())

    #beam = wolter3.trace_wolter3(beam, z0)
    beam = wolter3.trace_compound(beam)


    beam.plot_xz()

    beam.retrace(0.1)
    beam.plot_xz()

    plt.show()