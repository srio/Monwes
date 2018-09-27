from monwes.Beam import Beam
from monwes.Shape import BoundaryRectangle
import numpy as np
import matplotlib.pyplot as plt
from monwes.CompoundOpticalElement import CompoundOpticalElement
import h5py
from srxraylib.plot.gol import plot_scatter

do_plot = True
main = "__main__"


if main == "__main__":

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  example_montel_paraboloid")

    beam = Beam(25000)
    beam.set_circular_spot(1e-3)
    beam.set_flat_divergence(0.01, 0.01)
    beam.set_flat_divergence(1e-6, 1e-6)

    beam.plot_xz(0)

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


    montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=p, q=q, theta_z=theta, bound1=bound1, bound2=bound2, distance_of_the_screen=q)
    beam1, beam2, beam03 = montel.trace_montel(beam, name_file= 'polletto', print_footprint=0)



    f = h5py.File('polletto' + '.h5', 'r')
    n = np.ones(1)
    f['montel_good_rays/Number of rays'].read_direct(n)
    n = int(n[0])
    x1 = np.ones(n)
    y1 = np.ones(n)
    f['montel_good_rays/xoe1'].read_direct(x1)
    f['montel_good_rays/yoe1'].read_direct(y1)
    f.close()


    plot_scatter(y1*1e3, x1*1e3, title="footprint on oe1", xtitle='y[mm]', ytitle='x[mm]')