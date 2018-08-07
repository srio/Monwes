from Codes.Beam import Beam
from Codes.Shape import BoundaryRectangle
import unittest
import numpy as np
import matplotlib.pyplot as plt
from Codes.CompoundOpticalElement import CompoundOpticalElement


do_plot = True
main = "__main__"


if main == "__main__":

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_wolter_1_microscope")

    p = 13.4
    q =  0.300
    d = 0.1082882
    q1 = 0.67041707
    theta1 = 88.8*np.pi/180
    theta2 = 89.*np.pi/180

    wolter_jap = CompoundOpticalElement.wolter_for_japanese(p=p, q=q, d=d, q1=q1, theta1=theta1, theta2=theta2)


    beam = Beam()
    beam.set_gaussian_divergence(5*1e-5,0.00025)
    beam.set_rectangular_spot( xmax=200*1e-6, xmin=-200*1e-6, zmax=10*1e-6, zmin=-10*1e-6)

    beam.plot_xz(0)
    plt.title('wolter microscope')
    beam.plot_xpzp(0)


    wolter_jap.trace_wolter_japanese(beam)

    b2 = beam.y
    b3 = beam.z

    beam.y = b3
    beam.z = b2

    beam.plot_xz(0)
    beam.histogram()

    plt.show()