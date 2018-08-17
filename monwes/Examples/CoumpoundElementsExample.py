from monwes.Beam import Beam
from monwes.Shape import BoundaryRectangle
import numpy as np
import matplotlib.pyplot as plt
from numpy.testing import assert_almost_equal
from monwes.CompoundOpticalElement import CompoundOpticalElement

do_plot = True
main = "__main__"


def example_optimezed_wolter1_good_rays():
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_optimezed_wolter1_good_rays")

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

    assert_almost_equal(beam.x[indices], 0., 8)
    assert_almost_equal(beam.z[indices], 0., 8)

    beam.retrace(100.)
    beam.plot_good_xz()
    plt.title("optimezed_wolter1_good_rays")

    if do_plot:
        plt.show()


def example_wolter2_good_rays():
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_wolter2_good_rays")

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



def example_wolter3():
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


def example_kirk_patrick_baez():
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   example_kirk_patrick_baez")

    beam=Beam.initialize_as_person()
    beam.set_flat_divergence(1e-12, 1e-12)
    beam.x = beam.x*1e-3
    beam.z = beam.z*1e-3


    bound1 = BoundaryRectangle(xmax=2.5, xmin=-2.5, ymax=2.5, ymin=-2.5)
    bound2 = BoundaryRectangle(xmax=1., xmin=-1., ymax=1., ymin=-1.)

    kirk_patrick_baez = CompoundOpticalElement.initialize_as_kirkpatrick_baez(p=10., q=5., separation=4.,
                                                                              theta=89 * np.pi / 180, bound1=bound1,
                                                                              bound2=bound2)

    beam = kirk_patrick_baez.trace_compound(beam)

    beam.plot_good_xz(0)
    plt.title('Kirk Patrick Baez')

    indices = np.where(beam.flag > 0)

    beam.retrace(50.)

    beam.plot_good_xz()

    print(kirk_patrick_baez.info())

    print("Number of good rays: %f" % (beam.number_of_good_rays()))

    # beam.histogram()

    if do_plot:
        plt.show()


def example_wolter_1_microscope():

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


def example_montel_elliptical():
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  example_montel_elliptical")

    beam = Beam(25000)
    beam.set_flat_divergence(25*1e-6, 25*1e-6)
    beam.set_rectangular_spot(xmax=25*1e-6, xmin=-25*1e-6, zmax=5*1e-6, zmin=-5*1e-6)
    beam.set_gaussian_divergence(25*1e-4, 25*1e-4)



    beam.flag *= 0

    p = 5.
    q = 15.
    #theta = np.pi/2 - 0.15
    theta = 85. * np.pi / 180

    xmax = 0.
    xmin = -0.3
    ymax =  0.1
    ymin = -0.1
    zmax =  0.3
    zmin = 0.

    bound1 = BoundaryRectangle(xmax, xmin, ymax, ymin, zmax, zmin)
    bound2 = BoundaryRectangle(xmax, xmin, ymax, ymin, zmax, zmin)


    montel = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=p, q=q, theta=theta, bound1=bound1, bound2=bound2)
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

    plt.show()

def example_montel_paraboloid():
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




if main == "__main__":

    example_optimezed_wolter1_good_rays()
    example_wolter2_good_rays()
    example_wolter3()
    example_kirk_patrick_baez()
    example_wolter_1_microscope()
    example_montel_elliptical()
    example_montel_paraboloid()
