from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle

import numpy as np
import matplotlib.pyplot as plt

import Shadow

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





main = "__montel__ellipsoidal__"
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

    do_plot = True

    beam = Beam.initialize_from_shadow_beam(shadow_source())

    if do_plot:
        beam.plot_xz(0,title="source size")
        beam.plot_xpzp(0,title="source divergences")

    xmax = 0.
    xmin = -100.
    ymax = 0.3
    ymin = -0.4
    zmax = 100.
    zmin = 0.

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)

    montel_1 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.4, q=0.3, theta_z=theta,
                                                                     bound1=bound,
                                                                     bound2=bound, fp=0.4,
                                                                     fq=10000000.)

    #TODO force geometrical p and q
    beam_n, beam_m, beam_1 = montel_1.trace_montel(beam,do_plot_footprint=do_plot,name_file="tmp.h5")[2]


    if do_plot:
        beam_1.plot_xpzp(equal_axis=0,title='Plot of a the first montel ellipsoidal system')




    beam_2_04 = beam_1.duplicate()
    beam_2_14 = beam_1.duplicate()


    montel_2_04 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=0.4, theta_z=theta,
                                        bound1=bound, bound2=bound,
                                        fp=1000000000., fq=0.4)

    montel_2_14 = CompoundOpticalElement.initialize_as_montel_ellipsoid(p=0.3, q=1.471, theta_z=88.281 * np.pi / 180,
                                        bound1=bound, bound2=bound,
                                        fp=1000000000., fq=1.471)



    beam_n, beam_m, beam_2_04 = montel_2_04.trace_montel(beam_2_04,do_plot_footprint=do_plot)[2]
    #
    #
    if do_plot:
        beam_2_04.plot_xz(title='Final space plot of a montel ellipsoidal with parameter = 0.4')
        beam_2_04.plot_xpzp(title='Final velocity plot of a montel ellipsoidal with parameter = 0.4')

    # shadow_beam = beam_2_04.get_shadow_beam()


    beam_n, beam_m, beam_2_14 = montel_2_14.trace_montel(beam_2_14,do_plot_footprint=do_plot)[2]

    if do_plot:
        beam_2_14.plot_xz()
        plt.title('Final space plot of a montel ellipsoidal with parameter = 1.471')
        beam_2_14.plot_xpzp()
        plt.title('Final velocity plot of a montel ellipsoidal with parameter = 1.471')

    plt.show()



if main == "__montel__paraboloid__":

    mode = 1
    beam = beam()

    #beam = Beam(5000)
    #beam.set_flat_divergence(1e-4, 1e-5)
    #beam.set_rectangular_spot(xmax=1e-3, xmin=-1e-3, zmax=1e-4, zmin=-1e-4)


    #print("dx at the starting beam = %g\nIn the ideal cacse at the end of the system dx = %g" %(max(beam.x)-min(beam.x),max(beam.x)-min(beam.x)*1.471/0.4))
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

    montel_1 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.4, q=0.3, theta_z=theta, infinity_location='q')

    beam_n, beam_m, beam_1 = montel_1.trace_montel(beam, mode=mode)[2]

    beam_1.plot_xz()
    plt.title('Plot of a the first montel ellipsoidal system')

    beam_1.plot_xpzp()
    plt.title('Plot of a the divrgence first montel ellipsoidal system')

    beam_2_04 = beam_1.duplicate()
    beam_2_14 = beam_1.duplicate()


    montel_2_04 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.3, q=0.4, theta_z=theta, infinity_location='p')

    montel_2_14 = CompoundOpticalElement.initialize_as_montel_parabolic(p=0.3, q=1.471, theta_z=theta, infinity_location='p')


    beam_n, beam_m, beam_2_04 = montel_2_04.trace_montel(beam_2_04, mode=mode)[2]
    beam_n, beam_m, beam_2_14 = montel_2_14.trace_montel(beam_2_14, mode=mode)[2]


    #print("dx after 0.4 system = %g\n" %(max(beam_2_04.x)-min(beam_2_04.x)))
    #print("dx after 1.471 system = %g\n" %(max(beam_2_14.x)-min(beam_2_14.x)))


    beam_2_04.plot_xz(0)
    plt.title('Final plot of a montel paraboloid with parameter = 0.4')


    beam_2_14.plot_xz()
    plt.title('Final plot of a montel paraboloid with parameter = 1.471')


    plt.show()
