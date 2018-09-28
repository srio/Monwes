from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import time
import h5py
from srxraylib.plot.gol import plot_scatter

main = "__Moretti_1__"



def theta():

    p = 0.23e-3
    f = 0.23e-3/2
    xp = 0.0127046
    yp = 0.350885
    m = xp/p


    v = Vector(xp, yp-f, 0.)
    v.normalization()
    t = Vector(xp/p, -1, 0.)
    t.normalization()

    print((np.arccos(v.dot(t)))*180./np.pi)

    return np.arccos(v.dot(t))







if main == "__Moretti_1__":

    beam1 = Beam(1e5)
    beam1.set_gaussian_divergence(25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6, 25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6)

    xmax = 0.01
    xmin = -0.01
    ymax = 0.150
    ymin = -0.150
    zmax = 0.1
    zmin = -0.1

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
    dvx = [Beam(), Beam(), Beam(), Beam(), Beam()]
    bem = [Beam(), Beam(), Beam(), Beam(), Beam()]
    area = [Beam(), Beam(), Beam(), Beam(), Beam()]
    b = [Beam(), Beam(), Beam(), Beam(), Beam()]

    bins_list = [1, 2, 3, 4, 5]
    y_list = [1, 2, 3, 4, 5]


    Nn = 7


    for i in range (Nn):

        print("iteration %d" %i)

        beam = beam1.duplicate()
        alpha = round((-0.03 + 0.01 * i) * np.pi / 180 , 3)
        print(alpha)




        montel =  CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta_z=theta(), bound1=bound, bound2=bound, angle_of_mismatch=alpha, infinity_location='q')
        beam1, beam2, beam = montel.trace_montel(beam,print_footprint=0)[2]

        beam = b[i]



    b[0].plot_good_xpzp(0)