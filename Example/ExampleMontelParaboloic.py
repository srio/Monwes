from Codes.Beam import Beam
from Codes.Shape import BoundaryRectangle
import numpy as np
import matplotlib.pyplot as plt
from Codes.CompoundOpticalElement import CompoundOpticalElement

do_plot = True
main = "__main__"


if main == "__main__":

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  example_montel_paraboloid")

    beam = Beam(25000)
    beam.set_circular_spot(1e-3)
    beam.set_flat_divergence(0.01, 0.01)
    beam.set_flat_divergence(1e-6, 1e-6)


    beam.flag *= 0

    p = 5.
    q = 15.
    theta = 88.*np.pi/180

    xmax = 0.
    xmin = -0.4
    ymax =  0.4
    ymin = -0.4
    zmax =  0.4
    zmin = 0.

    bound1 = BoundaryRectangle(xmax, xmin, ymax, ymin, zmax, zmin)
    bound2 = BoundaryRectangle(xmax, xmin, ymax, ymin, zmax, zmin)


    montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=p, q=q, theta=theta, bound1=bound1, bound2=bound2, distance_of_the_screen=q)
    beam03 = montel.trace_montel(beam)

    print(beam03[2].N/25000)

    plt.figure()
    plt.plot(beam03[0].x, beam03[0].z, 'ro')
    plt.plot(beam03[1].x, beam03[1].z, 'bo')
    plt.plot(beam03[2].x, beam03[2].z, 'go')
    plt.xlabel('x axis')
    plt.ylabel('z axis')
    plt.axis('equal')

    beam03[2].plot_xz(0)

    print("No reflection = %d\nOne reflection = %d\nTwo reflection = %d" %(beam03[0].N, beam03[1].N, beam03[2].N))
    print("dx = %f" %(max(beam03[2].x)-min(beam03[2].x)))


    plt.show()
