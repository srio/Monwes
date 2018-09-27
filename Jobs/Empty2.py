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
import h5py
from srxraylib.plot.gol import plot_scatter
import Shadow


beam = Beam(5e4)
#beam.set_rectangular_spot(xmax=0.02e-2, xmin=-0.02e-2, zmax=0.001e-2, zmin=-0.001e-2)
beam.set_gaussian_divergence(4.6699999e-05, 0.00025499999)
beam.plot_xpzp(0)

op_axis = Beam(1)
op_axis.set_divergences_collimated()
beam = op_axis.merge(beam)

theta1 = 88.8 * np .pi / 180
theta2 = 89.  * np .pi / 180


d = 0.1082882
q = 0.3

ellipse = Optical_element.initialize_as_surface_conic_ellipsoid_from_focal_distances(p=13.4, q=0.,theta=theta1, fp=13.4, fq=0.67041707)
hyp = Optical_element.initialize_as_surface_conic_hyperboloid_from_focal_distances(p=0.562129, q=0.3, theta=theta2)

bound = BoundaryRectangle(xmax=0.75e-2, xmin=-0.75e-2, ymax=0.05, ymin=-0.05)
hyp.set_bound(bound)


theta_grazing1 = np.pi / 2 - theta1
theta_grazing2 = np.pi / 2 - theta2

print(beam.x[0], beam.y[0], beam.z[0])
print(d)
dy = d * np.cos(2 * theta_grazing1)
dz = d * np.sin(2 * theta_grazing1)


beam.y -= 13.4


ellipse.rotation_surface_conic(-theta_grazing1,'x')

[beam, t] = ellipse.intersection_with_optical_element(beam)
ellipse.output_direction_from_optical_element(beam)

print(beam.x[0], beam.y[0], beam.z[0])
print(beam.vx[0], beam.vy[0], beam.vz[0])
print(np.arctan(beam.vz[0]/beam.vy[0]) * 180. / np.pi)

#hyp.rotation_surface_conic(-np.pi , 'y')
hyp.rotation_surface_conic(-(2 * theta_grazing1 + theta_grazing2), 'x')
hyp.translation_surface_conic(dy, 'y')
hyp.translation_surface_conic(dz, 'z')
[beam, t] = hyp.intersection_with_optical_element(beam)
hyp.output_direction_from_optical_element(beam)


print(dy, dz)
print(beam.x[0], beam.y[0], beam.z[0])
print(beam.vx[0], beam.vy[0], beam.vz[0])
print(np.arctan(beam.vz[0]/beam.vy[0]) * 180. / np.pi)

beam.y -= dy
beam.z -= dz

beam.plot_zy(0)
plt.plot(beam.z[0], beam.y[0], 'k.')


position = Vector(beam.x,beam.y,beam.z)
velocity = Vector(beam.vx,beam.vy,beam.vz)
position.rotation(-(2 * theta_grazing1 + 2 *theta_grazing2),"x")
velocity.rotation(-(2 * theta_grazing1 + 2 *theta_grazing2),"x")
[beam.x,beam.y,beam.z] = [position.x,position.y,position.z]
[beam.vx,beam.vy,beam.vz] = [velocity.x,velocity.y,velocity.z]
print(beam.x[0], beam.y[0], beam.z[0])
print(beam.vx[0], beam.vy[0], beam.vz[0])

print(np.arctan(beam.vz[0]/beam.vy[0]) * 180. / np.pi)

beam.retrace(q)

beam.x *= 1e6
beam.z *= 1e6

beam.plot_good_xz(0, markersize=2.)
plt.title("final plot")

plt.show()