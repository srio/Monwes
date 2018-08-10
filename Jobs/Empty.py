from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
import os


def plot_montel(beam, equal_axis=1):
    plt.figure()
    plt.plot(beam[0].x, beam[0].z, 'r.')
    plt.plot(beam[1].x, beam[1].z, 'b.')
    plt.plot(beam[2].x, beam[2].z, 'g.')
    plt.xlabel('x axis')
    plt.ylabel('z axis')

    if equal_axis == 1:
        plt.axis('equal')

    beam[2].plot_xz(0)


    print("No reflection = %d\nOne reflection = %d\nTwo reflection = %d" %(beam[0].N, beam[1].N, beam[2].N))


beam = Beam()
#beam.set_circular_spot(1e-2)
beam.set_flat_divergence(1e-2, 1e-4)

beam.plot_xz(0)

index = np.where(beam.x != None)

xmax = 0.
xmin = -100.
ymax = 0.3
ymin = -0.4
zmax = 100.
zmin = 0.

bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

montel_1 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=5., q=25., theta=88.281 * np.pi / 180,
                                                                 bound1=bound,
                                                                 bound2=bound, distance_of_the_screen=25., fp=5.,
                                                                 fq=25.)

beam_3 = montel_1.trace_montel(beam)

plot_montel(beam_3)

montel_1.oe[0].output_frame_wolter(beam_3[2])


#
#unos = origin0
#
#dos = origin0*0.
#indices = np.where(unos==3)
#dos[indices] = dos[indices]*0.
#indices = np.where(unos==1)
#dos[indices] = origin1[0]
#indices = np.where(unos==2)
#dos[indices] = origin2[0]
#
#
#tros = origin0*0.
#indices = np.where(dos==3)
#tros[indices] = tros[indices]*0.
#indices = np.where(dos==1)
#tros[indices] = origin1[1]
#indices = np.where(dos==2)
#tros[indices] = origin2[1]
#
#
#x1 = np.ones(len(unos))*0.
#y1 = np.ones(len(unos))*0.
#z1 = np.ones(len(unos))*0.
#
#indices = np.where(unos==1)
#x1[indices] = beam_1[0].x
#y1[indices] = beam_1[0].y
#z1[indices] = beam_1[0].z
#indices = np.where(unos==2)
#x1[indices] = beam_2[0].x
#y1[indices] = beam_2[0].y
#z1[indices] = beam_2[0].z
#indices = np.where(unos==3)
#x1[indices] = beam_3[0].x
#y1[indices] = beam_3[0].y
#z1[indices] = beam_3[0].z
#indices = np.where(unos==0)
#x1[indices] = np.inf
#y1[indices] = np.inf
#z1[indices] = np.inf
#
#x2 = np.ones(len(unos))*0.
#y2 = np.ones(len(unos))*0.
#z2 = np.ones(len(unos))*0.
#
#indices = np.where(dos==1)
#x2[indices] = beam_1[1].x
#y2[indices] = beam_1[1].y
#z2[indices] = beam_1[1].z
#indices = np.where(dos==2)
#x2[indices] = beam_2[1].x
#y2[indices] = beam_2[1].y
#z2[indices] = beam_2[1].z
#indices = np.where(dos==3)
#x2[indices] = beam_3[1].x
#y2[indices] = beam_3[1].y
#z2[indices] = beam_3[1].z
#indices = np.where(dos==0)
#x2[indices] = np.inf
#y2[indices] = np.inf
#z2[indices] = np.inf
#
#indeces_x2_none = indices
#
#x3 = np.ones(len(unos))*0.
#y3 = np.ones(len(unos))*0.
#z3 = np.ones(len(unos))*0.
#
#indices = np.where(tros==1)
#x3[indices] = beam_1[2].x
#y3[indices] = beam_1[2].y
#z3[indices] = beam_1[2].z
#indices = np.where(tros==2)
#x3[indices] = beam_2[2].x
#y3[indices] = beam_2[2].y
#z3[indices] = beam_2[2].z
#indices = np.where(tros==3)
#x3[indices] = beam_3[2].x
#y3[indices] = beam_3[2].y
#z3[indices] = beam_3[2].z
#indices = np.where(tros==0)
#x3[indices] = np.inf
#y3[indices] = np.inf
#z3[indices] = np.inf
#
#
#index0 = np.array2string(index[0], threshold=np.inf)
#
#f = open("dati.txt", "w")
#
#
#
#f.write("Index     ons     dos     tros         x1          y1          z1        x2        y2         z2         x3         y3          z3\n")
#
#for i in range(0,len(x3)):
#
#    #print(x2[i])
#
#    if x1[i] == np.inf:
#        a1 = "     "
#    else:
#        a1 = ""
#
#
#    if x2[i] == np.inf:
#        a2 = "    "
#    else:
#        a2 = ""
#
#    if x3[i] == np.inf:
#        if x2[i] != np.inf:
#            a3 = "   "
#        else:
#            a3 = ""
#    else:
#        a3 = ""
#
#
#    if i < 10:
#        f.write(("%d           %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e2:
#        f.write(("%d          %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e3:
#        f.write(("%d         %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e4:
#        f.write(("%d        %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e5:
#        f.write(("%d       %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e6:
#        f.write(("%d      %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#    elif i < 1e7:
#        f.write(("%d     %d      %d       %d     %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s   %s%f%s\n" %(i,unos[i],dos[i],tros[i],a1,x1[i],a1,a1,y1[i],a1,a1,z1[i],a1,a2,x2[i],a2,a2,y2[i],a2,a2,z2[i],a2,a3,x3[i],a3,a3,y3[i],a3,a3,z3[i],a3)   ))
#
#
#f.close()
#
#a = "                  "
#
#print("%s 1" %a)
#
plt.show()