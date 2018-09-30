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

filename = "pollo"
title = "Rectangular spot of 1x1 um2 , gaussian divergence of FWHM 25 urad, displacement at (0, 0, 0)"





#theta = 88.9628519649*np.pi/180
#theta_grazing = np.pi/2 - theta
#print(88.9628519649)
#
#y1 = Vector(0., 1., 0.)
#y2 = Vector(0., 1., 0.)
#x = Vector(1., 0., 0.)
#z = Vector(0., 0., 1.)
#
#y1.rotation(theta_grazing,'x')
#print(np.arctan(y1.x / y1.y) * 180 / np.pi , np.arctan(y1.z / y1.y) * 180 / np.pi)
#
#y2.rotation(theta_grazing,'z')
#print(np.arctan(y2.x / y2.y) * 180 / np.pi + 90. , np.arctan(y2.z / y2.y) * 180 / np.pi)
#
#y3 = Vector(y2.x, y1.y, y1.z)
#y3.normalization()
#print(np.arctan(y3.x / y3.y) * 180 / np.pi + 90. , np.arctan(y3.z / y3.y) * 180 / np.pi - 90)
#print(y3.info())
#
#y4 = Vector(y3.x, y3.y, 0.)
#y4.normalization()
#y5 = Vector(0., y3.y, y3.z)
#y5.normalization()
#
#y4.rotation(-theta_grazing, 'z')
#y5.rotation(-theta_grazing, 'x')
#print("Eddaje")
#print(y4.info())
#print(y5.info())
#
#
#
#
#
#
#print('\n')
#
#vf = 250.
#
#beam = Beam()
#beam.set_circular_spot(1e-3)
##beam.set_rectangular_spot(xmin=-0.5e-9, xmax=0.5e-9, zmin=-0.5e-9, zmax=0.5e-9)
#beam.set_gaussian_divergence(vf* 1e-6 / (2 * np.sqrt(2 * np.log(2))), vf * 1e-6 / (2 * np.sqrt(2 * np.log(2))))
##beam.set_gaussian_spot(1e-6 / (2 * np.sqrt(2 * np.log(2))), 1.e-6 / (2 * np.sqrt(2 * np.log(2))))
#beam.set_divergences_collimated()
#
##op_axis = Beam(1)
##op_axis.set_point(0., 0., 0.)
##op_axis.set_divergences_collimated()
##beam = op_axis.merge(beam)
#
#beam.plot_xz(0)
#beam.plot_xpzp(0)
#
#
#print(beam.vx, beam.vz)
#
#xmax = 0.0
#xmin = -0.01
#ymax = 0.300
#ymin = -0.300
#zmax = 0.01
#zmin = 0.0
#
#bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
#
#
#
#montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=88.9628519649*np.pi/180, bound1=bound, bound2=bound, infinity_location='p')
#par = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.351, q=1., theta=theta, infinity_location='q')
#
#
#
#
#
#
#beam = montel.trace_montel(beam)[2]
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(np.arctan(beam.vx[0] / beam.vy[0]) * 180 / np.pi + 90. , 90.- np.arctan(beam.vz[0] / beam.vy[0]) * 180 / np.pi)
#
#shadow_beam = Shadow.Beam(beam.N)
#
#shadow_beam.rays[:, 0] = beam.x
#shadow_beam.rays[:, 1] = beam.y
#shadow_beam.rays[:, 2] = beam.z
#shadow_beam.rays[:, 3] = beam.vx
#shadow_beam.rays[:, 4] = beam.vy
#shadow_beam.rays[:, 5] = beam.vz
#
#shadow_beam.rays[:, 6] = np.zeros(beam.N) + 1.0
#shadow_beam.rays[:, 7] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 8] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 9] = np.zeros(beam.N) + 1.0  # beam_2_04.flag
#shadow_beam.rays[:, 10] = np.zeros(beam.N) + 1e-8
#shadow_beam.rays[:, 11] = np.arange(beam.N) + 1
#
#shadow_beam.rays[:, 12] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 13] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 14] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 15] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 16] = np.zeros(beam.N) + 0.0
#shadow_beam.rays[:, 17] = np.zeros(beam.N) + 0.0
#
#shadow_beam.write("beam_ollo.sha")
#
#
#beam.plot_xz(0)
#beam.plot_xpzp(0)
#


theta = 88.*np.pi/180
theta_grazing = np.pi/2 - theta

print("Start Vector")

y = Vector(0., 1., 0.)
y.z = - np.tan(theta_grazing) / np.sqrt(1+np.tan(theta_grazing)**2)
y.x = + np.tan(theta_grazing) / np.sqrt(1+np.tan(theta_grazing)**2)
y.y =  np.sqrt(1 - y.x**2 - y.z**2)

y.x = y.x
y.y = y.y
y.z = y.z
print(y.x, y.y, y.z)

print(np.arctan(y.x / np.sqrt(y.y ** 2 + y.z**2)) * 180 / np.pi - 90.,   90. + np.arctan(y.z / np.sqrt(y.y ** 2 + y.x ** 2)) * 180 / np.pi)

alpha = -np.arctan(y.z/y.y)
y.rotation(alpha, 'x')
gamma = np.arctan(y.x/y.y)
y.rotation(gamma, 'z')
print(y.x, y.y, y.z)

y = Vector(0., 1., 0.)
y.rotation(-gamma, 'z')
y.rotation(-alpha, 'x')
print(y.x, y.y, y.z)

print("End Vector")

#y = Vector(0., 1., 0.)
#y.z = - np.tan(theta_grazing) / np.sqrt(2 + (np.tan(theta_grazing))**2)
#y.x =  np.tan(theta_grazing) / np.sqrt(2 + (np.tan(theta_grazing))**2)
#y.y =  1 / np.sqrt(2 + np.tan(theta_grazing)**2)
#y.normalization()
#print(y.x, y.y, y.z)
#print(np.arctan(y.x/y.y)*180/np.pi, np.arctan(y.z/y.y)*180/np.pi)
#
#alpha = -np.arctan(y.z/y.y)
#y.rotation(alpha, 'x')
#gamma = np.arctan(y.x/y.y)
#y.rotation(gamma, 'z')
#
#print(y.x, y.y, y.z)



beam = Beam()
beam.set_rectangular_spot(xmax=0.5e-3, xmin=-0.5e-3, zmax=0.5e-3, zmin=-0.5e-3)
beam.set_flat_divergence(500*1e-6,800*1e-6)
beam.set_divergences_collimated()

op_axis = Beam(1)
op_axis.set_point(0., 0., 0.)
op_axis.set_divergences_collimated()
beam = op_axis.merge(beam)


xmax = 0.0
xmin = -0.01
ymax = 0.300
ymin = -0.300
zmax = 0.01
zmin = 0.0
bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta=theta, bound1=bound, bound2=bound, infinity_location='p')
par = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(0.351,1.,theta,0.0,'p',None,1)



#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-gamma, 'z')
#velocity.rotation(-alpha, 'x')
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#
#
#print(beam.vx, beam.vy, beam.vz)
#
#
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(-gamma, 'z')
#position.rotation(-alpha, 'x')
#beam.x = position.x
#beam.y = position.y
#beam.z = position.z
#
#
#beam.x -= beam.vx[0] * 0.351
#beam.y -= beam.vy[0] * 0.351
#beam.z -= beam.vz[0] * 0.351
#
#
#print(beam.x, beam.y, beam.z)
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(np.arctan(beam.vx[0] / np.sqrt(beam.vy[0] ** 2 + beam.vz[0] ** 2)) * 180 / np.pi , np.arctan(beam.vz[0] / np.sqrt(beam.vy[0] ** 2 + beam.vx[0] ** 2)) * 180 / np.pi)
#
#beam = montel.trace_montel(beam)[2]
#
##[beam, t] = par.intersection_with_optical_element(beam)
##par.output_direction_from_optical_element(beam)
#
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(np.arctan(beam.vx[0] / np.sqrt(beam.vy[0] ** 2 + beam.vz[0] ** 2)) * 180 / np.pi , np.arctan(beam.vz[0] / np.sqrt(beam.vy[0] ** 2 + beam.vx[0] ** 2)) * 180 / np.pi)
#
#
#
#
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-alpha, 'x')
#velocity.rotation(-gamma, 'z')
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#
#beam.plot_xpzp(0)
#
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(-alpha, 'x')
#position.rotation(-gamma, 'z')
#beam.x = position.x
#beam.y = position.y
#beam.z = position.z
#print(beam.x[0], beam.y[0], beam.z[0])
#
#beam.retrace(1.)
#
#beam.plot_xz(0)
#beam.plot_xpzp(0)
#
#plt.show()







#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-gamma, 'z')
#velocity.rotation(-alpha, 'x')
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#
#
#
#
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(-gamma, 'z')
#position.rotation(-alpha, 'x')
#beam.x = position.x
#beam.y = position.y
#beam.z = position.z
#
#
#beam.x -= beam.vx[0] * 0.351
#beam.y -= beam.vy[0] * 0.351
#beam.z -= beam.vz[0] * 0.351
#
#
#print(np.arctan(beam.vx[0]/beam.vy[0])*180/np.pi, np.arctan(beam.vz[0]/beam.vy[0])*180/np.pi)
#
##[beam, t] = par.intersection_with_optical_element(beam)
##par.output_direction_from_optical_element(beam)
#
#print(np.arctan(beam.vx[0]/beam.vy[0])*180/np.pi, np.arctan(beam.vz[0]/beam.vy[0])*180/np.pi)
#
#beam = montel.trace_montel(beam)[2]
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(np.arctan(beam.vx[0] / np.sqrt(beam.vy[0] ** 2 + beam.vz[0] ** 2)) * 180 / np.pi - 90.,90. + np.arctan(beam.vz[0] / np.sqrt(beam.vy[0] ** 2 + beam.vx[0] ** 2)) * 180 / np.pi)
#print(np.arctan(beam.vx[0]/beam.vy[0])*180/np.pi, np.arctan(beam.vz[0]/beam.vy[0]))
#
#
#
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-alpha, 'x')
#velocity.rotation(-gamma, 'z')
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(-gamma, 'z')
#position.rotation(-alpha, 'x')
#beam.x = position.x
#beam.y = position.y
#beam.z = position.z
#print(beam.x[0], beam.y[0], beam.z[0])
#
#beam.retrace(1.)
#
#beam.plot_xz(0)
#beam.plot_xpzp(0)


def tras_rot(beam, p):
    velocity = Vector(beam.vx, beam.vy, beam.vz)
    velocity.rotation(-gamma, 'z')
    velocity.rotation(-alpha, 'x')
    velocity.normalization()
    beam.vx = velocity.x
    beam.vy = velocity.y
    beam.vz = velocity.z

    position = Vector(beam.x, beam.y, beam.z)
    position.rotation(-gamma, 'z')
    position.rotation(-alpha, 'x')
    beam.x = position.x
    beam.y = position.y
    beam.z = position.z

    beam.x -= beam.vx[0] * p
    beam.y -= beam.vy[0] * p
    beam.z -= beam.vz[0] * p


def out(beam, q):
    print("\nOut")
    velocity = Vector(beam.vx, beam.vy, beam.vz)
    velocity.rotation(-alpha, 'x')
    print(velocity.x[0], velocity.y[0], velocity.z[0])
    velocity.rotation(-gamma, 'z')
    print(velocity.x[0], velocity.y[0], velocity.z[0])
    beam.vx = velocity.x
    beam.vy = velocity.y
    beam.vz = velocity.z

    position = Vector(beam.x, beam.y, beam.z)
    position.rotation(-alpha, 'x')
    position.rotation(-gamma, 'z')
    beam.x = position.x
    beam.y = position.y
    beam.z = position.z

    beam.retrace(q)





beam = Beam(25000)
beam.set_flat_divergence(0.1e-6, 0.1e-6)
#beam.set_rectangular_spot(xmax=0.00015, xmin=-0.00015, zmax=0.00008, zmin=-0.008)
#beam.set_divergences_collimated()
#beam.plot_xpzp(0)
op_axis = Beam(1)
op_axis.set_point(0., 0., 0.)
op_axis.set_divergences_collimated()
beam = op_axis.merge(beam)


xmax = 0.
xmin = -100.
ymax = 0.3
ymin = -0.4
zmax = 100.
zmin = 0.

montel_1 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.4, q=0.3, theta=theta, bound1=bound, bound2=bound, infinity_location='q')
par = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3,  theta=theta, cylindrical=1, infinity_location='q')
par1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3, theta=theta, cylindrical=1, infinity_location='q')

#beam.plot_xz(0)
#beam.plot_xpzp(0)
#tras_rot(beam, 0.351)
#beam.vx -= beam.vx[0]
#beam.vy -= beam.vy[0]
#beam.vz -= beam.vz[0]
#beam.plot_xz(0)
#beam.plot_xpzp(0)
#print(beam.vx, beam.vy, beam.vz)


#par.rotation_surface_conic(-theta_grazing, 'x')
#beam.y -= 0.4
#
#
#[beam, t] = par.intersection_with_optical_element(beam)
#par.output_direction_from_optical_element(beam)
#print(np.arctan(beam.vz[0]/beam.vy[0])*180/np.pi)
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-2*theta_grazing, 'x')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#
#beam.plot_xpzp(0)



#beam.y -= 0.4
#par.rotation_surface_conic(-np.pi/2, 'y')
##par.rotation_surface_conic(theta_grazing, 'z')
#[beam, t] = par.intersection_with_optical_element(beam)
#print(beam.x[0], beam.y[0], beam.z[0])
#par.output_direction_from_optical_element(beam)
#print(np.arctan(beam.vx[0]/beam.vy[0])*180/np.pi)
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(2*theta_grazing, 'z')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#beam.plot_xpzp(0)





def time_comparison(oe, beam1, elements):
    origin = np.ones(beam1.N)
    tf = 1e35 * np.ones(beam1.N)
    tau = np.ones(3) * 6

    for i in range(0, len(elements)):
        beam = beam1.duplicate()

        [beam, t] = oe[i].intersection_with_optical_element(beam)
        # t = abs(t)
        tau[i] = np.mean(t)

        indices = np.where(beam.flag < 0)

        t[indices] = 1e30

        tf = np.minimum(t, tf)
        indices = np.where(t == tf)
        origin[indices] = elements[i]

    return origin


par2 = par1.duplicate()
par2.rotation_surface_conic(-theta_grazing, 'x')
par1.rotation_surface_conic(-theta_grazing, 'x')


ccc = np.array([0., 0., 0., 0., 0., 0., 0., 0., 1., 0.])
par2 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
par2.set_parameters(0.4, 0.3, 0., 0., "Surface conical mirror")
par2.rotation_surface_conic(-theta_grazing, 'x')



#par1.rotation_surface_conic(2*theta_grazing, 'x')
par2.rotation_surface_conic(np.pi/2, 'y')
#par2.rotation_surface_conic(theta_grazing, 'z')

#ccc = np.array([0., 0., 0., 0., 0., 0., 0., 1., 0., -10.])
#screen = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
#screen.set_parameters(0.4, 10., 0., 0., "Surface conical mirror")
#
#beam.y -= 0.4
#
#origin = time_comparison([par1, par2, screen], beam, [1,2,3])
#
#indices = np.where(origin == 1)
#beam1 = beam.part_of_beam(indices)
#indices = np.where(origin == 2)
#beam2 = beam.part_of_beam(indices)
#indices = np.where(origin == 3)
#beam3 = beam.part_of_beam(indices)
#
#print(beam1.N, beam2.N, beam3.N)
#
#if beam3.N != 0:
#    [beam3, t] = screen.intersection_with_optical_element(beam3)
#
#beam1_list = [beam1.duplicate(), Beam(), Beam()]
#beam2_list = [beam2.duplicate(), Beam(), Beam()]
#beam3_list = [beam3.duplicate(), Beam(), Beam()]
#
#origin1 = [1, 2]
#origin2 = [1, 2]
#
#
#print(beam1_list[0].vx[2], beam1_list[0].vy[2], beam1_list[0].vz[2])
#
#[beam1_list[1], t] = par2.intersection_with_optical_element(beam2_list[0])
#
#
#[beam1_list[1], t] = par1.intersection_with_optical_element(beam1_list[1])
#
#
##velocity = Vector(beam1_list[1].vx, beam1_list[1].vy, beam1_list[1].vz)
##velocity.rotation(-theta_grazing, 'x')
##print(velocity.z)
##theta_plot = -np.arctan(velocity.z /  np.sqrt(velocity.y**2 + velocity.x**2)) * 180 / np.pi
###theta_plot = -np.arctan(velocity.z /  np.sqrt(velocity.y**2)) * 180 / np.pi
##y = beam1_list[1].y
#
#par1.output_direction_from_optical_element(beam1_list[1])
#
#beam2_list[2] = beam1_list[1]
#
#
##print(beam2_list[2].vx[2], beam2_list[2].vy[2], beam2_list[2].vz[2])
##
##[beam2_list[2], t] = par1.intersection_with_optical_element(beam1_list[1])
##
##velocity = Vector(beam2_list[2].vx, beam2_list[2].vy, beam2_list[2].vz)
##velocity.rotation(-theta_grazing, 'x')
##print(velocity.z)
##theta_plot = -np.arctan(velocity.z / np.sqrt(velocity.y**2 + velocity.x**2)) * 180 / np.pi
###theta_plot = -np.arctan(velocity.z /  np.sqrt(velocity.y**2)) * 180 / np.pi
##y = beam2_list[2].y
##
##
##par1.output_direction_from_optical_element(beam2_list[2])
##
##print(beam2_list[2].vx[2], beam2_list[2].vy[2], beam2_list[2].vz[2])
##
#
#velocity = Vector(beam2_list[2].vx, beam2_list[2].vy, beam2_list[2].vz)
##velocity.rotation(np.arctan(velocity.x[0]/velocity.y[0]), 'z')
##print(velocity.x[0], velocity.y[0], velocity.z[0])
##velocity.rotation(-np.arctan(velocity.z[0]/velocity.y[0]), 'x')
##print(velocity.x[0], velocity.y[0], velocity.z[0])
##beam2_list[2].vx = velocity.x
##beam2_list[2].vy = velocity.y
##beam2_list[2].vz = velocity.z
#
#
#beam2_list[2].plot_xpzp(0)
#plt.title("final")
#
#
#print(velocity.modulus())
#
##plot_scatter(y, theta_plot, title="2 -> 1 (bad one) v2 ", xtitle="y", ytitle="theta")
#
#
#plt.show()




def tras_rot2(beam,p,theta):
    theta_grazing = np.pi / 2 -theta
    y = Vector(0., 1., 0.)
    y.rotation(-theta_grazing, 'x')
    x = Vector(1., 0., 0.)
    xp = x.vector_product(y)
    vrot = y.rodrigues_formula(xp,-theta_grazing)
    vrot.normalization()

    position = Vector(beam.x, beam.y, beam.z)
    velocity = Vector(beam.vx, beam.vy, beam.vz)
    position.rotation(-theta_grazing, 'x')
    velocity.rotation(-theta_grazing, 'x')
    #position = position.rodrigues_formula(xp, -theta_grazing)
    #velocity = velocity.rodrigues_formula(xp, -theta_grazing)
    velocity.normalization()
    beam.vx = velocity.x
    beam.vy = velocity.y
    beam.vz = velocity.z

    vector_point = Vector(0, p, 0)
    vector_point.rotation(-theta_grazing, "x")
    print(vector_point.info())
    #vector_point = vector_point.rodrigues_formula(xp, -theta_grazing)
    vector_point.normalization()

    beam.x = position.x - vector_point.x * p
    beam.y = position.y - vector_point.y * p
    beam.z = position.z - vector_point.z * p



#beam = Beam(25e3)
#beam.set_flat_divergence(10e-6, 10e-6)
#
#
#par = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3, theta=88*np.pi/180, cylindrical=1, infinity_location='q')
#par.rotation_surface_conic(-theta, 'z')
#
#ccc = np.array([0., 0., 0., 0., 0., 0., 1., 0., 0., 0.])
#oe2 = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
#oe2.set_parameters(0.4, 0.3, 0., 0., "Surface conical mirror")
#
#
##par.rotation_surface_conic(np.pi/2, 'y')
#
#ccc = np.array([0., 0., 0., 0., 0., 0., 0., 1., 0., -10])
#screen = Optical_element.initialize_as_surface_conic_from_coefficients(ccc)
#screen.set_parameters(0.4, 0.3, 0., 0., "Surface conical mirror")
#
#
#tras_rot(beam, 0.4)
#print(beam.x[0], beam.y[0], beam.z[0])
#
#
#origin = time_comparison([par, oe2, screen], beam, [1,2,3])
#indices = np.where(origin == 1)
#beam1 = beam.part_of_beam(indices)
#indices = np.where(origin == 2)
#beam2 = beam.part_of_beam(indices)
#indices = np.where(origin == 3)
#beam3 = beam.part_of_beam(indices)
#
#
#print(beam1.N, beam2.N, beam3.N)
#
#if beam3.N != 0:
#    [beam3, t] = screen.intersection_with_optical_element(beam3)
#
#beam1_list = [beam1.duplicate(), Beam(), Beam()]
#beam2_list = [beam2.duplicate(), Beam(), Beam()]
#beam3_list = [beam3.duplicate(), Beam(), Beam()]
#
#[beam2_list[1], t] = par.intersection_with_optical_element(beam2_list[0])
#par.output_direction_from_optical_element(beam2_list[1])
#
#velocity = Vector(beam2_list[1].vx, beam2_list[1].vy, beam2_list[1].vz)
#velocity.rotation(np.arctan(velocity.x[0]/velocity.y[0]), 'z')
#velocity.rotation(-np.arctan(velocity.z[0]/velocity.y[0]), 'x')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#beam2_list[1].vx = velocity.x
#beam2_list[1].vy = velocity.y
#beam2_list[1].vz = velocity.z


#[beam2_list[2], t] = par2.intersection_with_optical_element(beam2_list[1])
#print(beam2_list[2].x, beam2_list[2].y, beam2_list[2].z)
#oe1.output_direction_from_optical_element(beam2_list[2])

#beam2_list[1].plot_xpzp(0)


#beam = Beam(25000)
#beam.set_flat_divergence(200e-6, 200e-6)
#beam.set_rectangular_spot(xmax=0.00015, xmin=-0.00015, zmax=0.00008, zmin=-0.008)
#beam.set_divergences_collimated()
#
#beam.plot_xpzp(0)
#beam.plot_xz(0)
#
#theta = 88. * np.pi / 180
#
#xmax = 0. + 1e-9
#xmin = -100.
#ymax = 0.3
#ymin = -0.4
#zmax = 100.
#zmin = 0. - 1e-9
#
#bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
#
#montel = CompoundOpticalElement.initialize_as_new_montel_paraboloid(p=0.4, q=0.3, theta=theta, infinity_location='p')
#
#beam = montel.trace_new_montel(beam)[2]
#
#beam.retrace(0.3)
#
##velocity = Vector(beam.vx, beam.vy, beam.vz)
##alpha = np.arctan(velocity.x[0]/velocity.y[0])
##velocity.rotation(alpha, 'z')
##beta = np.arctan(velocity.z[0]/velocity.y[0])
##velocity.rotation(-beta, 'x')
##beam.vx = velocity.x
##beam.vy = velocity.y
##beam.vz = velocity.z
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(beam.vz)
#
#beam.plot_xpzp(0)
#beam.plot_xz(0)
#
#plt.show()




#beam = Beam(25000)
#beam.set_flat_divergence(0.1e-1, 0.1e-1)
##beam.set_rectangular_spot(xmax=0.00015, xmin=-0.00015, zmax=0.00008, zmin=-0.008)
##beam.set_divergences_collimated()
#beam.plot_xpzp(0)
#op_axis = Beam(1)
#op_axis.set_point(0., 0., 0.)
#op_axis.set_divergences_collimated()
#beam = op_axis.merge(beam)
#
#theta = 80. * np.pi / 180
#theta_grazing = np.pi / 2 - theta
#
#xmax = 0.
#xmin = -100.
#ymax = 0.3
#ymin = -0.4
#zmax = 100.
#zmin = 0.
#
#montel_1 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.4, q=0.3, theta=theta, bound1=bound, bound2=bound, infinity_location='q')
#par = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3,  theta=theta, cylindrical=1, infinity_location='q')
#par1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3, theta=theta, cylindrical=1, infinity_location='q')
#par2 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.4, q=0.3, theta=theta, cylindrical=1, infinity_location='q')
#
#
#
##vx = np.tan(theta_grazing) / np.sqrt(1 + np.tan(theta_grazing)**2)
##y = Vector(vx, np.sqrt(1-2*vx**2), -vx)
##
##print(y.info())
##
##
##alpha = np.arctan(y.x/y.y)
##y.rotation(alpha, 'z')
##beta = np.arctan(y.z/y.y)
##y.rotation(-beta, 'x')
##
##print(y.info())
##
##
##velocity = Vector(beam.vx, beam.vy, beam.vz)
##velocity.rotation(-alpha, 'z')
##velocity.rotation(beta, 'x')
##beam.vx = velocity.x
##beam.vy = velocity.y
##beam.vz = velocity.z
##
##
##position = Vector(beam.x, beam.y, beam.z)
##position.rotation(-alpha, 'z')
##position.rotation(beta, 'x')
##beam.x = position.x - velocity.x[0]
##beam.y = position.y - velocity.y[0]
##beam.z = position.z - velocity.z[0]
#
#y = Vector(0., 1., 0.)
#y.rotation(-theta_grazing, 'z')
#vz = np.sin(theta_grazing)
#vx = y.x * np.sqrt(1. - np.sin(theta_grazing)**2) / np.sqrt(y.y**2 + y.x**2)
#vy = np.sqrt(1. - vx**2 - vz**2)
#y2 = Vector(vx, vy, vz)
#
#
#x = Vector(1., 0., 0.)
#x.rotation(-theta_grazing, 'z')
#y = y.rodrigues_formula(x, theta_grazing)
#print(y.info())
#
#
#y2 = Vector(0., 1., 0.)
#y2 = y2.rodrigues_formula(x, theta_grazing)
#alpha = np.arctan(y2.z / y2.y)
#beta = np.arctan(y2.y / y2.x)
#print(y2.info())
#
#y3 = Vector(0., 1., 0.)
#y3.rotation(alpha, 'x')
#print(y3.info())
#
#print(alpha*180/np.pi , beta*180/np.pi)
#
#
#print(np.arctan(y.z / np.sqrt(y.x**2 + y.y **2)) * 180. / np.pi)
#
#
##velocity = Vector(beam.vx, beam.vy, beam.vz)
###velocity.rotation(-theta_grazing, 'x')
###print(np.arctan(velocity.z[0]/ np.sqrt(velocity.x[0]**2 + velocity.y[0]**2)) * 180. / np.pi)
##velocity.rotation(-theta_grazing, 'z')
##velocity = velocity.rodrigues_formula(x, -theta_grazing)
##beam.vx = velocity.x
##beam.vy = velocity.y
##beam.vz = velocity.z
##position = Vector(beam.x, beam.y, beam.z)
###position.rotation(-theta_grazing, 'x')
##position.rotation(-theta_grazing, 'z')
##position = position.rodrigues_formula(x, -theta_grazing)
##beam.x = position.x - velocity.x[0] * 0.4
##beam.y = position.y - velocity.y[0] * 0.4
##beam.z = position.z - velocity.z[0] * 0.4
##
##print(beam.vx[0], beam.vy[0], beam.vz[0])
##print(beam.x[0], beam.y[0], beam.z[0])
##
##par.rotation_surface_conic(theta_grazing, 'z')
##
##
##print(np.arctan(beam.vz[0]/np.sqrt(beam.vx[0]**2+beam.vy[0]**2)) + 90.)
##print(np.arctan(beam.vx[0]/np.sqrt(beam.vx[0]**2+beam.vy[0]**2)) + 90.)
##
##print(90. - np.arctan(beam.vz[0]/beam.vy[0]))
##print(90. - np.arctan(beam.vx[0]/beam.vy[0]))
##
##[beam, t] = par.intersection_with_optical_element(beam)
##par.output_direction_from_optical_element(beam)
##
##
##print(beam.vx[2], beam.vy[2], beam.vz[2])
##
##
##
##velocity = Vector(beam.vx, beam.vy, beam.vz)
##velocity = velocity.rodrigues_formula(x, -theta_grazing)
##velocity.rotation(theta_grazing, 'z')
###velocity.rotation(-theta_grazing, 'x')
###alpha = np.arctan(velocity.x[0]/velocity.y[0])
##print(velocity.x[0], velocity.y[0], velocity.z[0])
##beam.vx = velocity.x
##beam.vy = velocity.y
##beam.vz = velocity.z
###
###
###
###print(beam.vx[2], beam.vy[2], beam.vz[2])
###
##beam.plot_xpzp(0)
##plt.show()
#
#
#
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(-theta_grazing, 'z')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#velocity = velocity.rodrigues_formula(x, -theta_grazing)
##velocity.rotation(-theta_grazing, 'x')
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(-theta_grazing, 'z')
#position = position.rodrigues_formula(x, -theta_grazing)
##position.rotation(-theta_grazing, 'x')
#beam.x = position.x - velocity.x[0] * 0.4
#beam.y = position.y - velocity.y[0] * 0.4
#beam.z = position.z - velocity.z[0] * 0.4
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(beam.x[0], beam.y[0], beam.z[0])
#
#
#par.rotation_surface_conic(np.pi/2, 'y')
#par.rotation_surface_conic(-alpha, 'x')
#par.rotation_surface_conic(-beta, 'y')
#
#[beam, t] = par.intersection_with_optical_element(beam)
#par.output_direction_from_optical_element(beam)
#
#
#print(beam.vx[2], beam.vy[2], beam.vz[2])
#
#print(np.arctan(beam.vz[0]/np.sqrt(beam.vx[0]**2+beam.vy[0]**2)))
#print(np.arctan(beam.vx[0]/np.sqrt(beam.vx[0]**2+beam.vy[0]**2)))
#
#print(np.arctan(beam.vz[0]/beam.vy[0]))
#print(np.arctan(beam.vx[0]/beam.vy[0]))
##
##velocity = Vector(beam.vx, beam.vy, beam.vz)
##velocity = velocity.rodrigues_formula(z, -theta_grazing)
##print(velocity.x[0], velocity.y[0], velocity.z[0])
##velocity.rotation(0*theta_grazing, 'x')
##velocity.rotation(-theta_grazing, 'z')
##alpha = np.arctan(velocity.x[0]/velocity.y[0])
##print(velocity.x[0], velocity.y[0], velocity.z[0])
##beam.vx = velocity.x
##beam.vy = velocity.y
##beam.vz = velocity.z
##
##
###print(beam.vx[2], beam.vy[2], beam.vz[2])
##
#
#
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(np.arctan(velocity.x[0]/velocity.y[0]), 'z')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#velocity.rotation(-np.arctan(velocity.z[0]/velocity.y[0]), 'x')
#print(velocity.x[0], velocity.y[0], velocity.z[0])
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#
#beam.plot_xpzp(0)
#
#
#
#plt.show()


#
#
#beam = Beam(25000)
#beam.set_gaussian_divergence(10.6*1e-6, 10.6*1e-6)
#beam.set_rectangular_spot(xmax=0.5e-4, xmin=-0.5e-4, zmax=0.5e-4, zmin=-0.5e-4)
#beam.set_circular_spot(1e-4)
#
#
#beam.plot_xpzp(0)
#op_axis = Beam(1)
#op_axis.set_point(0., 0., 0.)
#op_axis.set_divergences_collimated()
#beam = op_axis.merge(beam)
#
#theta = 88. * np.pi / 180
#theta_grazing = np.pi / 2 - theta
#
#xmax = 0.
#xmin = -100.
#ymax = 0.3
#ymin = -0.4
#zmax = 100.
#zmin = 0.
#
#
#
#
#bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
#montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1, theta=theta, bound1=bound, bound2=bound, infinity_location='q', angle_of_mismatch=0.0)
#
#vector = Vector(0., 1., 0.)
#vector.rotation(-theta_grazing, 'x')
#x = Vector(1., 0., 0.)
#n = x.vector_product(vector)
#
#
#vrot = vector.rodrigues_formula(n, -theta_grazing)
#vrot.normalization()
#print(vrot.info())
#
#print(np.arctan(vrot.x/vrot.y)*180/np.pi, np.arctan(vrot.z/vrot.y)*180/np.pi)
#print(np.arctan(vrot.x/np.sqrt(vrot.y**2+vrot.z**2))*180/np.pi, np.arctan(vrot.z/np.sqrt(vrot.x**2+vrot.y**2))*180/np.pi)
#
#
#vector = Vector(0., 1., 0.)
#vector.rotation(-theta_grazing, 'z')
#x = Vector(0., 0., 1.)
#n = x.vector_product(vector)
#vrot = vector.rodrigues_formula(n, -theta_grazing)
#print(vrot.info())
#
#print(np.arctan(vrot.x/vrot.y)*180/np.pi, np.arctan(vrot.z/vrot.y)*180/np.pi)
#print(np.arctan(vrot.x/np.sqrt(vrot.y**2+vrot.z**2))*180/np.pi, np.arctan(vrot.z/np.sqrt(vrot.x**2+vrot.y**2))*180/np.pi)
#






beam = Beam(25000)
beam.set_gaussian_divergence(10.6*1e-6, 10.6*1e-6)
#beam.set_rectangular_spot(xmax=0.5e-4, xmin=-0.5e-4, zmax=0.5e-4, zmin=-0.5e-4)
#beam.set_circular_spot(1e-4)


beam.plot_xpzp(0)
op_axis = Beam(1)
op_axis.set_point(0., 0., 0.)
op_axis.set_divergences_collimated()
beam = op_axis.merge(beam)

theta = 88. * np.pi / 180
theta_grazing = np.pi / 2 - theta

xmax = 0.
xmin = -100.
ymax = 0.3
ymin = -0.4
zmax = 100.
zmin = 0.




bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
montel = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1, theta=theta, bound1=bound, bound2=bound, infinity_location='q', angle_of_mismatch=0.0)


#vx = np.tan(theta_grazing)/np.sqrt(1+np.tan(theta_grazing)**2)
#
#vector = Vector(vx, np.sqrt(1-2*vx**2), -vx)
#print(vector.info())
#alpha = np.arctan(vector.x/vector.y)
#vector.rotation(alpha, 'z')
#beta = np.arctan(vector.z/vector.y)
#vector.rotation(-beta, 'x')
#print(vector.info())
#
#
#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(beta, 'x')
#velocity.rotation(-alpha, 'z')
#
#position = Vector(beam.x, beam.y, beam.z)
#position.rotation(beta, 'x')
#position.rotation(-alpha, 'z')
#
#
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z
#
#
#beam.x = position.x - velocity.x[0]*0.351
#beam.y = position.y - velocity.y[0]*0.351
#beam.z = position.z - velocity.z[0]*0.351
#
#print(beam.vx[0], beam.vy[0], beam.vz[0])
#print(beam.x[0], beam.y[0], beam.z[0])
#
#print(np.arctan(beam.vx[0]/beam.vy[0])*180/np.pi, np.arctan(beam.vz[0]/beam.vy[0])*180/np.pi)
#print(np.arctan(beam.vx[0]/np.sqrt(beam.vy[0]**2+beam.vz[0]**2))*180/np.pi, np.arctan(beam.vz[0]/np.sqrt(beam.vx[0]**2+beam.vy[0]**2))*180/np.pi)


beam = montel.trace_montel(beam)[2]


#velocity = Vector(beam.vx, beam.vy, beam.vz)
#velocity.rotation(np.arctan(np.mean(velocity.x)/np.mean(velocity.y)), 'z')
#print(np.mean(velocity.x), np.mean(velocity.y), np.mean(velocity.z))
#velocity.rotation(-np.arctan(np.mean(velocity.z)/np.mean(velocity.y)), 'x')
#print(np.mean(velocity.x), np.mean(velocity.y), np.mean(velocity.z))
#
#
#
#beam.vx = velocity.x
#beam.vy = velocity.y
#beam.vz = velocity.z

beam.plot_xpzp(0)

plt.show()


