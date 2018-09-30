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

main = "__Moretti_2__"



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

    beam1 = Beam(1e6)
    #beam1.set_gaussian_spot(1 / (2 * np.sqrt(2 * np.log(2))) * 1e-6, 1 / (2 * np.sqrt(2 * np.log(2))) * 1e-6)
    #beam1.set_rectangular_spot(xmax=0.5e-6, xmin=-0.5e-6, zmax=0.5e-6, zmin=-0.5e-6)
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

    bins_list = [1, 2, 3, 4, 5]
    y_list = [1, 2, 3, 4, 5]

    date = "dati/Moretti/" + "M.Moretti  " + time.strftime("%d-%m-%y at %H:%M:%S") + ".h5"

    f = h5py.File(date, 'w')
    f.close()

    Nn = 7


    for i in range (Nn):


        beam = beam1.duplicate()
        alpha = round(-0.03 + 0.01 * i , 3)
        print(alpha)



        f = h5py.File(date, 'a')
        f1 = f.create_group(str(alpha))

        montel =  CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta(), bound1=bound, bound2=bound, angle_of_mismatch=alpha)
        beam = montel.trace_montel(beam,f1)


        f1["Number of rays"] = beam[2].N
        f1["x"] = beam[2].x
        f1["y"] = beam[2].y
        f1["z"] = beam[2].z
        f1["vx"] = beam[2].vx
        f1["vy"] = beam[2].vy
        f1["vz"] = beam[2].vz

        f.close()


    #beam = beam1.duplicate()
    #alpha = 9.600000
    #print(alpha)

    #f = h5py.File(date, 'a')
    #f1 = f.create_group(str(alpha))

    #montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta(), bound1=bound,
    #                                                               bound2=bound, angle_of_mismatch=alpha)
    #beam = montel.trace_montel(beam, f1)

    #f1["Number of rays"] = beam[2].N
    #f1["x"] = beam[2].x
    #f1["y"] = beam[2].y
    #f1["z"] = beam[2].z
    #f1["vx"] = beam[2].vx
    #f1["vy"] = beam[2].vy
    #f1["vz"] = beam[2].vz

    #f.close()


    plt.show()


if main == '__Moretti_2__':

    dim = 0.5e-6

    for ii in range (0,1):

        vf = 25

        print("The value of divergence fwhm is %f" %vf)

        beam1 = Beam()
        #beam1.set_gaussian_spot(dim / (2 * np.sqrt(2 * np.log(2))) * 1e-6, dim / (2 * np.sqrt(2 * np.log(2))) * 1e-6)
        #beam1.set_gaussian_divergence(vf / (2 * np.sqrt(2 * np.log(2))) * 1e-6, 25 / (vf * np.sqrt(2 * np.log(2))) * 1e-6)

        #beam1 = Beam()
        #beam1.set_rectangular_spot(xmax=0.5e-6, xmin=-0.5e-6, zmax=0.5e-6, zmin=-0.5e-6)
        #beam1.set_circular_spot(1e-6/np.pi)
        beam1.set_gaussian_divergence(vf


                                      * 1e-6, vf / (2 * np.sqrt(2 * np.log(2))) * 1e-6)



        alpha_initial = -0.02
        alpha_step = 0.0005

        print("alpha initial = %f, alpha step = %f" %(alpha_initial, alpha_step))


        xmax = 0.01
        xmin = -0.01
        ymax = 0.150
        ymin = -0.150
        zmax = 0.01
        zmin = -0.01

        bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
        dvx = [Beam(), Beam(), Beam(), Beam(), Beam()]
        bem = [Beam(), Beam(), Beam(), Beam(), Beam()]
        area = [Beam(), Beam(), Beam(), Beam(), Beam()]

        bins_list = [1, 2, 3, 4, 5]
        y_list = [1, 2, 3, 4, 5]

        date = "dati/Moretti/" + "M.Moretti  " + time.strftime("%d-%m-%y at %H:%M:%S") + ".h5"

        f = h5py.File(date, 'w')
        f.close()

        Nn = 80


        ####### zero angle  #################################################################################################

        f = h5py.File(date, 'a')

        beam = beam1.duplicate()
        print(beam.N)
        f1 = f.create_group(str(0.00))
        montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta(), bound1=bound,
                                                                       bound2=bound, angle_of_mismatch=0.0)
        beam = montel.trace_montel(beam, f1)
        f1["Number of rays"] = beam[2].N
        f1["x"] = beam[2].x
        f1["y"] = beam[2].y
        f1["z"] = beam[2].z
        f1["vx"] = beam[2].vx
        f1["vy"] = beam[2].vy
        f1["vz"] = beam[2].vz
        f.close()

        ####################################################################################################################


        for i in range(Nn):
            beam = beam1.duplicate()
            alpha = round(alpha_initial+ alpha_step * i, 5)

            if abs(alpha) > 1e-13:
                print(alpha)

                f = h5py.File(date, 'a')
                f1 = f.create_group(str(alpha))

                montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta(), bound1=bound,
                                                                               bound2=bound, angle_of_mismatch=alpha)


                beam = montel.trace_montel(beam, f1)

                f1["Number of rays"] = beam[2].N
                f1["x"] = beam[2].x
                f1["y"] = beam[2].y
                f1["z"] = beam[2].z
                f1["vx"] = beam[2].vx
                f1["vy"] = beam[2].vy
                f1["vz"] = beam[2].vz

                f.close()


