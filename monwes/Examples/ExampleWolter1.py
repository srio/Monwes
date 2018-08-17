from monwes.Beam import Beam
import numpy as np
import matplotlib.pyplot as plt
from monwes.CompoundOpticalElement import CompoundOpticalElement

do_plot = True
main = "__main__"


if main == "__main__":

    p = 100.
    beam1 = Beam.initialize_as_person()
    beam1.x *= 50.
    beam1.z *= 50.
    beam1.set_point(p, 0., p)
    op_ax = Beam(1)
    op_ax.set_point(p, 0., p)
    beam = op_ax.merge(beam1)
    beam.set_divergences_collimated()
    beam.plot_xz()

    p = 1e12
    R = 100.
    theta = 1e-3 * np.pi / 180

    wolter1 = CompoundOpticalElement.initialiaze_as_wolter_1_with_two_parameters(p1=p, R=R, theta=theta)

    beam = wolter1.trace_good_rays(beam)
    beam.plot_good_xz()

    indices = np.where(beam.flag >= 0)

    beam.retrace(100.)
    beam.plot_good_xz()
    plt.title("optimezed_wolter1_good_rays")

    if do_plot:
        plt.show()
