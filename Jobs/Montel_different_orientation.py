from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime as dt
import time


#
#  main== 'focalize' it mean that the script focalize a beam with a circular source of 0.1 mm geometry and collimating divergence
#  main== 'collimate' it mean that the script collimate a beam with pointwise source geometry and gaussian divergence with sigma of 25e-6
#  mode 0 correspond to have vx=vz=tan(theta_grazing)/(1+tan(theta_grazing**2))
#  mode 1 correspond to have vx=vz such that tan(theta_grazing)=vx/vy
#  mode 2 do a first rotation on x and the second in a plane orthogonal to x vector y'
#  mode 3 do a first rotation on z and the second in a plane orthogonal to z vector y'
#  beam_sourrce1() is for the "focalize" beam
#  beam_source2() is for the "collimate" beam


main = 'collimate'
mode = 1

__name__ = '__main__'

def beam_source1():

    beam = Beam(25000)
    beam.set_circular_spot(r=1e-4)
    #beam.set_rectangular_spot(xmax=1e-4, xmin=-1e-4, zmax=1e-4, zmin=1e-4)
    #beam.set_gaussian_spot(dx=1e-4, dz=1e-4)
    #beam.set_flat_divergence(dx=25e-6, dz=25e-6)
    #beam.set_gaussian_divergence(25e-6, 25e-6)
    beam.set_divergences_collimated()
    xmax = 0.
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    return beam, bound


def beam_source2():

    beam = Beam(25000)
    #beam.set_circular_spot(r=1e-4)
    #beam.set_rectangular_spot(xmax=1e-4, xmin=-1e-4, zmax=1e-4, zmin=1e-4)
    #beam.set_gaussian_spot(dx=1e-4, dz=1e-4)
    #beam.set_flat_divergence(dx=25e-6, dz=25e-6)
    beam.set_gaussian_divergence(25e-6, 25e-6)
    #beam.set_divergences_collimated()
    xmax = 0.
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    return beam, bound


if __name__ == '__main__':

    if main == 'focalize':
        beam, bound = beam_source1()

        p = 0.4
        q = 1.
        theta = 88. * np.pi / 180

        beam.plot_xz(0)
        plt.title('Initial source dimension')

        montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=p, q=q, theta=theta, infinity_location='p', bound1=bound, bound2=bound)

        beam = montel.trace_montel(beam, mode=mode)[2]

        beam.plot_xz()
        title = 'Beam dimension after montel of mode = ' + str(mode)
        plt.title(title)



    if main == 'collimate':
        beam, bound = beam_source2()

        p = 0.4
        q = 1.
        theta = 88. * np.pi / 180


        beam.plot_xpzp(0)
        plt.title('Initial source divergence')

        montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=p, q=q, theta=theta,
                                                                       infinity_location='q', bound1=bound,
                                                                       bound2=bound)

        beam = montel.trace_montel(beam, mode=mode)[2]



        beam.plot_xpzp(0)
        title = 'Beam divergence after montel of mode = ' + str(mode)
        plt.title(title)


    plt.show()