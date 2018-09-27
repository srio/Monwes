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

p = 100.
beam1 = Beam.initialize_as_person()
beam1.x *= 50.
beam1.z *= 50.
beam1.set_point(p, 0., p)
beam1.plot_xz()
op_ax = Beam(1)
op_ax.set_point(p, 0., p)
beam = op_ax.merge(beam1)
beam.set_divergences_collimated()
beam.plot_xz()

p = 1e12
R = 100.
theta = 89. * np.pi / 180

wolter1 = CompoundOpticalElement.initialiaze_as_wolter_1_with_two_parameters(p1=p, R=R, theta=theta)

beam = wolter1.trace_good_rays(beam)
beam.plot_good_xz()

indices = np.where(beam.flag >= 0)

beam.retrace(100.)
beam.plot_good_xz(0)
plt.title("optimezed_wolter1_good_rays")
beam.plot_yx(0)

print(np.arctan(beam.vz[0] / beam.vy[0]) * 180. / np.pi)
print(wolter1.type)


print(beam.flag)
indices = np.where(beam.flag>=0)
v = Vector(np.mean(beam.vx[indices]), -np.mean(beam.vy[indices]), np.mean(beam.vz[indices]))


velocity = Vector(beam.vx, beam.vy, beam.vz)
position = Vector(beam.x, beam.y, beam.z)
alpha = np.arctan(v.x/v.y)
v.rotation(alpha, 'z')
velocity.rotation(alpha, 'z')
position.rotation(alpha, 'z')
alpha = np.arctan(-v.z/v.y)
velocity.rotation(alpha)
position.rotation(alpha)
v.rotation(alpha, 'x')
velocity.rotation(np.pi)

print(v.info())

beam.x = position.x
beam.y = position.y
beam.z = position.z


beam.vx = velocity.x
beam.vy = velocity.y
beam.vz = velocity.z

print(beam.vx, beam.vy, beam.vz)
beam.retrace(44)
beam.plot_xz(0)
plt.plot("fujihsa")

plt.show()