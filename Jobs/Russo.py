from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
#import Shadow

##################  Instruction #######################
# The different possible choises are:
#
# main == "__ideal_lenses__"             system composed by two ideal lensen in series
# main == "__KirkPatrickBaez__"          system composed by two kirkPatrikBaez system with parabolic mirrors
# main == "__paraboloid__"               system composed by two paraboloid mirror with no cylindricity
# main == "__montel__ellipsoidal__"      system composed by two montel system composed by ellipsoidal mirrors
# main == "__montel__paraboloid__"       system composed by two montel system composed by ellipsoidal paraboloid
#
# The source utilized is that imported by the oasys python script equal to the one of the "Echo Spectrometer"





main = "__montel__paraboloid__"
theta = 88.281*np.pi/180

def shadow_source():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 1
    oe0.FSOUR = 1
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.005
    oe0.HDIV2 = 0.005
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 0
    oe0.NPOINT = 25000
    oe0.PH1 = 9137.65
    oe0.VDIV1 = 0.00075
    oe0.VDIV2 = 0.00075
    oe0.WXSOU =  0.0003
    oe0.WZSOU =  0.00015

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")
    return beam

def beam():

    shadow_beam = shadow_source()

    beam = Beam(25000)
    beam.initialize_from_arrays(
        shadow_beam.getshonecol(1),
        shadow_beam.getshonecol(2),
        shadow_beam.getshonecol(3),
        shadow_beam.getshonecol(4),
        shadow_beam.getshonecol(5),
        shadow_beam.getshonecol(6),
        shadow_beam.getshonecol(10),
    )

    beam.flag *= 0.


    beam.plot_xz(0)
    plt.title('starting beam')

    return beam


if main == "__ideal_lenses__":

    beam = beam()

    ideal_lens_1 = Optical_element.initialiaze_as_ideal_lens(p=0.4, q=0.3, fx=0.4, fz=0.4)
    ideal_lens_2_04 = Optical_element.initialiaze_as_ideal_lens(p=0.3, q=0.4, fx=0.4, fz=0.4)
    ideal_lens_2_1471 = Optical_element.initialiaze_as_ideal_lens(p=0.3, q=1.471, fx = 1.471, fz = 1.471)


    beam = ideal_lens_1.trace_optical_element(beam)

    beam_04 = beam.duplicate()
    beam_1471 = beam.duplicate()


    beam_04 = ideal_lens_2_04.trace_optical_element(beam_04)        ########### case of 0.4
    beam_1471 = ideal_lens_2_1471.trace_optical_element(beam_1471)  ########### case of 1.471

    beam_04.plot_xz()
    plt.title('Final plot of an ideal lense with parameter = 0.4')
    beam_1471.plot_xz()
    plt.title('Final plot of an ideal lense with parameter = 1.471')

    plt.show()



if main == "__KirkPatrickBaez__":

    beam = beam()

    comp_oe1 = CompoundOpticalElement.initialize_as_kirkpatrick_baez_parabolic(p=0.4, q=0.3, separation=0.2, theta=88.281*np.pi/180, infinity_location='q')

    beam = comp_oe1.trace_compound(beam)


    beam_2_04 = beam.duplicate()
    beam_2_14 = beam.duplicate()

    comp_oe2_04 = CompoundOpticalElement.initialize_as_kirkpatrick_baez_parabolic(p=0.3, q=0.4, separation=0.2, theta=theta, infinity_location='p')
    comp_oe2_14 = CompoundOpticalElement.initialize_as_kirkpatrick_baez_parabolic(p=0.3, q=1.471, separation=0.2, theta=theta, infinity_location='p')

    beam_2_04 = comp_oe2_04.trace_compound(beam_2_04)
    beam_2_14 = comp_oe2_14.trace_compound(beam_2_14)


    beam_2_04.plot_xz()
    plt.title('Final plot of a KirkPatrickBaez system with parameter = 0.4')
    beam_2_14.plot_xz()
    plt.title('Final plot of a KirkPatrickBaez system with parameter = 1.471')


    plt.show()



if main == "__paraboloid__":

    beam = beam()

    oe1 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.3, q=0.4, theta=theta, infinity_location='q')

    beam = oe1.trace_optical_element(beam)


    beam_2_04 = beam.duplicate()
    beam_2_14 = beam.duplicate()

    oe2_04 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.3, q=0.4, theta=theta, infinity_location='p')
    oe2_14 = Optical_element.initialize_as_surface_conic_paraboloid_from_focal_distances(p=0.3, q=1.471, theta=theta, infinity_location='p')

    beam_2_04 = oe2_04.trace_optical_element(beam_2_04)
    beam_2_14 = oe2_14.trace_optical_element(beam_2_14)


    beam_2_04.plot_xz()
    plt.title('Final plot of a paraboloid mirrors with parameter = 0.4')
    beam_2_14.plot_xz()
    plt.title('Final plot of a paraboloid mirrors with parameter = 1.471')


    plt.show()



if main == "__montel__ellipsoidal__":

    beam = beam()

    beam.plot_xpzp(0)
    beam.plot_xz(0)



    xmax = 0.
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_1 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.4, q=0.3, theta=theta,
                                                                     bound1=bound,
                                                                     bound2=bound, fp=0.4,
                                                                     fq=10000000.)

    beam_1 = montel_1.trace_montel(beam)[2]

    beam_1.plot_xpzp(0)
    plt.plot(beam_1.x[0], beam_1.z[0], 'ko')
    plt.title('Plot of a the first montel ellipsoidal system')

    beam_1.plot_xpzp()



    beam_2_04 = beam_1.duplicate()
    beam_2_14 = beam_1.duplicate()


    montel_2_04 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=0.4, theta=theta,
                                                                     bound1=bound, bound2=bound, fp=1000000000., fq=0.4)

    montel_2_14 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=1.471, theta=88.281 * np.pi / 180,
                                                                     bound1=bound, bound2=bound,
                                                                     fp=1000000000., fq=1.471)


    beam_2_04 = montel_2_04.trace_montel(beam_2_04)[2]


    beam_2_04.plot_xz()
    plt.plot(beam_2_04.x[0], beam_2_04.z[0], 'ko')
    plt.title('Final space plot of a montel ellipsoidal with parameter = 0.4')
    beam_2_04.plot_xpzp()
    plt.plot(beam_2_04.vx[0]*1e6, beam_2_04.vz[0]*1e6, 'ko')
    plt.title('Final velocity plot of a montel ellipsoidal with parameter = 0.4')

    #shadow
    shadow_beam = Shadow.Beam(beam_2_04.N)

    shadow_beam.rays[:,0] = beam_2_04.x
    shadow_beam.rays[:,1] = beam_2_04.y
    shadow_beam.rays[:,2] = beam_2_04.z
    shadow_beam.rays[:,3] = beam_2_04.vx
    shadow_beam.rays[:,4] = beam_2_04.vy
    shadow_beam.rays[:,5] = beam_2_04.vz

    shadow_beam.rays[:, 6]  = np.zeros(beam_2_04.N) + 1.0
    shadow_beam.rays[:, 7]  = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 8]  = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 9]  = np.zeros(beam_2_04.N) + 1.0 # beam_2_04.flag
    shadow_beam.rays[:, 10] = np.zeros(beam_2_04.N) + 1e-8
    shadow_beam.rays[:, 11] = np.arange(beam_2_04.N) + 1

    shadow_beam.rays[:, 12] = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 13] = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 14] = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 15] = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 16] = np.zeros(beam_2_04.N) + 0.0
    shadow_beam.rays[:, 17] = np.zeros(beam_2_04.N) + 0.0



    shadow_beam.write("beam_2_04.sha")


    shadow_beam2 = Shadow.Beam()
    shadow_beam2.load("beam_2_04.sha")
    #Shadow.ShadowTools.plotxy(shadow_beam2, 1, 3)




    beam_2_14 = montel_2_14.trace_montel(beam_2_14)[2]


    beam_2_14.plot_xz()
    plt.plot(beam_2_14.x[0], beam_2_14.z[0], 'ko')
    plt.title('Final space plot of a montel ellipsoidal with parameter = 1.471')
    beam_2_14.plot_xpzp()
    plt.plot(beam_2_14.vx[0]*1e6, beam_2_14.vz[0]*1e6, 'ko')
    plt.title('Final velocity plot of a montel ellipsoidal with parameter = 1.471')

    plt.show()



if main == "__montel__paraboloid__":

    #beam = beam()

    beam = Beam(5000)
    beam.set_flat_divergence(1e-4, 1e-5)
    beam.set_rectangular_spot(xmax=1e-3, xmin=-1e-3, zmax=1e-4, zmin=-1e-4)


    print("dx at the starting beam = %g\nIn the ideal cacse at the end of the system dx = %g" %(max(beam.x)-min(beam.x),max(beam.x)-min(beam.x)*1.471/0.4))
    beam.plot_xpzp(0)
    plt.title("Initial divergence")
    beam.plot_xz(0)
    plt.title("Initial dimension")

    xmax = 0.
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_1 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.4, q=0.3, theta=theta, bound1=bound, bound2=bound, infinity_location='q')

    beam_1 = montel_1.trace_montel(beam)[2]

    beam_1.plot_xz()
    plt.plot(beam_1.x[0], beam_1.z[0], 'ko')
    plt.title('Plot of a the first montel ellipsoidal system')

    beam_1.plot_xpzp()
    plt.title('Plot of a the divrgence first montel ellipsoidal system')

    beam_2_04 = beam_1.duplicate()
    beam_2_14 = beam_1.duplicate()


    montel_2_04 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.3, q=0.4, theta=theta, bound1=bound, bound2=bound, infinity_location='p')

    montel_2_14 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.3, q=1.471, theta=theta, bound1=bound, bound2=bound, infinity_location='p')


    beam_2_04 = montel_2_04.trace_montel(beam_2_04)[2]
    beam_2_14 = montel_2_14.trace_montel(beam_2_14)[2]


    print("dx after 0.4 system = %g\n" %(max(beam_2_04.x)-min(beam_2_04.x)))
    print("dx after 1.471 system = %g\n" %(max(beam_2_14.x)-min(beam_2_14.x)))


    beam_2_04.plot_xz(0)
    plt.plot(beam_2_04.x[0], beam_2_04.z[0], 'ko')
    plt.title('Final plot of a montel paraboloid with parameter = 0.4')

    beam_2_04.retrace(1.)

    beam_2_04.plot_xz(0)
    plt.title('Final plot of a montel paraboloid with parameter = 0.4')

    beam_2_14.plot_xz()
    plt.plot(beam_2_14.x[0], beam_2_14.z[0], 'ko')
    plt.title('Final plot of a montel paraboloid with parameter = 1.471')


    plt.show()
