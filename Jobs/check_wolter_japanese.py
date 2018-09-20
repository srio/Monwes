from monwes.Beam import Beam
from monwes.Shape import BoundaryRectangle
import numpy as np
import matplotlib.pyplot as plt
from monwes.CompoundOpticalElement import CompoundOpticalElement

do_plot = True
main = "__main__"


if main == "__main__":


    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_wolter_1_microscope")

    p = 13.4
    q =  0.300
    d = 1.                            #0.1082882
    q1 = 5.                           #0.67041707
    theta1 = 88.8*np.pi/180
    theta2 = 88.8*np.pi/180

    wolter_jap = CompoundOpticalElement.wolter_1_for_microscope(p=p, q=q, d=d, q1=q1, theta1=theta1, theta2=theta2)


    beam = Beam()
    beam.set_gaussian_divergence(5*1e-5,0.00025)
    beam.set_rectangular_spot( xmax=200*1e-6, xmin=-200*1e-6, zmax=10*1e-6, zmin=-10*1e-6)
    #beam.set_divergences_collimated()

    beam.plot_xz(0)
    plt.title('wolter microscope')
    beam.plot_xpzp(0)

    op_axis = Beam(1)
    op_axis.set_point(0., 0., 0.)
    op_axis.set_divergences_collimated()

    beam = op_axis.merge(beam)


    beam = wolter_jap.trace_compound(beam)

    #b2 = beam.y
    #b3 = beam.z

    #beam.y = b3
    #beam.z = b2

    #beam.x *= 1e6
    #beam.z *= 1e6

    #beam.plot_xz(0)
    #beam.histogram()

    print("\n\n Outside the trace")
    print("Position")
    print(beam.x, beam.y, beam.z)
    print("Velocity")
    print(beam.vx, beam.vy, beam.vz)

    t = - beam.y/beam.vy
    beam.x += beam.vx * t
    beam.y += beam.vy * t
    beam.z += beam.vz * t

    beam.plot_xz(0)




    print("\n\n At the focus of hyperbola")
    print("Position")
    print(beam.x, beam.y, beam.z)
    print("Velocity")
    print(beam.vx, beam.vy, beam.vz)

    #wolter_jap.oe[0].new_output_frame_montel(beam)


    #Nn = 10000
    #dx = np.ones(Nn)
    #t = np.ones(Nn)

    #for i in range (Nn):

    #    dt = 0.4/Nn
    #    if i == 0:
    #        t[i] = 0
    #    else:
    #        t[i] = t[i-1] + dt
    #    beam.retrace(dt)

    #    dx[i] = max(beam.x) - min(beam.x)

    #plt.figure()
    #plt.plot(t,dx)
    #plt.ylabel('dx')
    #plt.xlabel('distance after the second mirror intersection')

    #beam.retrace(-0.4)

    #m = np.where(dx==min(dx))
    #beam.retrace(t[m])

    #beam.retrace(1.)

    #beam.plot_xz(0)
    #beam.plot_xpzp(0)

    plt.show()

