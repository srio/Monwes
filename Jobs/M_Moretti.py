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
import os

main = "__Moretti__"


def theta():

    p = 0.23e-3
    f = 0.23e-3/2
    xp = 0.012908941406853667
    yp = 0.3622625396643068
    m = xp/p


    v = Vector(xp, yp-f, 0.)
    v.normalization()
    t = Vector(1., m, 0.)
    t.normalization()

    print((np.pi/2 - np.arccos(v.dot(t)))*180./np.pi)

    return np.pi/2 - np.arccos(v.dot(t))





if main == "__Moretti__":

    beam1 = Beam()
    beam1.set_rectangular_spot(xmax=0.5e-6, xmin=-0.5e-6, zmax=0.5e-6, zmin=-0.5e-6)
    beam1.set_gaussian_divergence(25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6, 25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6)

    xmax = 0.
    xmin = -0.4
    ymax = 0.150
    ymin = -0.150
    zmax = 0.4
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
    dvx = [Beam(), Beam(), Beam(), Beam(), Beam(), Beam()]
    bem = [Beam(), Beam(), Beam(), Beam(), Beam(), Beam()]
    area = [Beam(), Beam(), Beam(), Beam(), Beam(), Beam()]

    bins_list = [1, 2, 3, 4, 5, 6]
    y_list = [1, 2, 3, 4, 5, 6]


    for i in range (6):

        beam = beam1.duplicate()
        alpha = -0.03 + 0.01 * i
        montel =  CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta(), bound1=bound, bound2=bound, angle_of_mismatch=alpha)
        beam = montel.trace_montel(beam)
        dvx[i] = beam[2].vx * 1e6
        bem[i] = beam[2]

        #bvx = dvx[i] / np.sum(dvx[i])
        area[i] = beam[2].N

        (mu, sigma) = norm.fit(dvx[i])
        n, bins , patches = plt.hist(dvx[i], 50)
        plt.title(alpha*-1)
        y = mlab.normpdf(bins, mu, sigma)

        bins_list[i] = bins
        y_list[i] = y

        plt.close()



    plt.figure()
    plt.plot(bins_list[0], y_list[0]*area[0], 'r')
    plt.plot(bins_list[1], y_list[1]*area[1], 'b')
    plt.plot(bins_list[2], y_list[2]*area[2], 'g')
    plt.plot(bins_list[3], y_list[3]*area[3], 'k')
    plt.plot(bins_list[4], y_list[4]*area[4], 'y')
    plt.plot(bins_list[5], y_list[5]*area[5], 'c')
    plt.legend([0.03,0.02,0.01,0.00,-0.01,-0.02])

    plt.figure()
    plt.hist(dvx[0], histtype='step', normed=0, color='r')
    plt.hist(dvx[1], histtype='step', normed=0, color='b')
    plt.hist(dvx[2], histtype='step', normed=0, color='g')
    plt.hist(dvx[3], histtype='step', normed=0, color='k')
    plt.hist(dvx[4], histtype='step', normed=0, color='y')
    plt.hist(dvx[5], histtype='step', normed=0, color='c')
    plt.legend([0.03,0.02,0.01,0.00,-0.01,-0.02])



    plt.figure()
    plt.plot(bem[0].vx, bem[0].vz, 'r.')
    plt.title("vx/vz plot of 0.03 degree")

    plt.figure()
    plt.plot(bem[1].vx, bem[1].vz, 'b.')
    plt.title("vx/vz plot of 0.02 degree")

    plt.figure()
    plt.plot(bem[2].vx, bem[2].vz, 'g.')
    plt.title("vx/vz plot of 0.01 degree")

    plt.figure()
    plt.plot(bem[3].vx, bem[3].vz, 'k.')
    plt.title("vx/vz plot of 0.00 degree")

    plt.figure()
    plt.plot(bem[4].vx, bem[4].vz, 'y.')
    plt.title("vx/vz plot of -0.01 degree")

    plt.figure()
    plt.plot(bem[5].vx, bem[5].vz, 'c.')
    plt.title("vx/vz plot of -0.02 degree")




    plt.figure()
    plt.plot(bem[0].x, bem[0].z, 'r.')
    plt.title("x/z plot of 0.03 degree")

    plt.figure()
    plt.plot(bem[1].x, bem[1].z, 'b.')
    plt.title("x/z plot of 0.02 degree")

    plt.figure()
    plt.plot(bem[2].x, bem[2].z, 'g.')
    plt.title("x/z plot of 0.01 degree")

    plt.figure()
    plt.plot(bem[3].x, bem[3].z, 'k.')
    plt.title("x/z plot of 0.00 degree")

    plt.figure()
    plt.plot(bem[4].x, bem[4].z, 'y.')
    plt.title("x/z plot of -0.01 degree")

    plt.figure()
    plt.plot(bem[5].x, bem[5].z, 'c.')
    plt.title("x/z plot of -0.02 degree")



    plt.show()
