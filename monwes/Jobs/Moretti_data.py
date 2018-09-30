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

save = False
show = False
x_shift = True

main = '__Moretti_5__'
dim = 'big'

filename = 'M.Moretti  07-09-18 at 16:46:53'
filename = 'dati/Moretti/' + filename + '.h5'

alpha_initial = - 0.02
alpha_step = 0.0005


title = "(-10um,0,0)"

def import_beam(string = None):

    n = np.ones(1)
    f[string + '/Number of rays'].read_direct(n)


    beam = Beam(int(n[0]))

    if string is not None:

        f[string + '/x'].read_direct(beam.x)
        f[string + '/y'].read_direct(beam.y)
        f[string + '/z'].read_direct(beam.z)
        f[string + '/vx'].read_direct(beam.vx)
        f[string + '/vy'].read_direct(beam.vy)
        f[string + '/vz'].read_direct(beam.vz)



    return beam

def import_footprint(string = None):

    n = np.ones(1)
    f[string + '/Number of rays'].read_direct(n)


    xoe1 = np.ones(int(n[0]))
    yoe1 = np.ones(int(n[0]))
    zoe2 = np.ones(int(n[0]))
    yoe2 = np.ones(int(n[0]))
    on_do = np.ones(int(n[0]))

    f[string + '/montel_good_rays/xoe1'].read_direct(xoe1)
    f[string + '/montel_good_rays/yoe1'].read_direct(yoe1)
    f[string + '/montel_good_rays/zoe2'].read_direct(zoe2)
    f[string + '/montel_good_rays/yoe2'].read_direct(yoe2)
    f[string + '/montel_good_rays/on_do'].read_direct(on_do)

    footprint = [xoe1, yoe1, zoe2, yoe2, on_do]

    return footprint



def plot_footprint(footprint, angle):


    xoe1 =  footprint[0]
    yoe1 =  footprint[1]
    zoe2 =  footprint[2]
    yoe2 =  footprint[3]
    on_do = footprint[4]

    plt.figure()
    indices = np.where(on_do==1)
    plt.plot(yoe1[indices]*1e6,xoe1[indices]*1e6, color='r', marker='.',markersize=0.4, linestyle='None')
    indices = np.where(on_do==2)
    plt.plot( yoe1[indices]*1e6,xoe1[indices]*1e6, color='b', marker='.',markersize=0.4, linestyle='None')
    plt.title('footprint on oe1 of angle ∆=%.2f°' %angle)
    plt.xlabel('y[um]')
    plt.ylabel('x[um]')

    print(on_do)

    plt.figure()
    indices = np.where(on_do==1)
    plt.plot(yoe2[indices]*1e6,zoe2[indices]*1e6, color='r', marker='.',markersize=0.4, linestyle='None')
    indices = np.where(on_do==2)
    plt.plot( yoe2[indices]*1e6,zoe2[indices]*1e6, color='b', marker='.',markersize=0.4, linestyle='None')
    plt.title('footprint on oe2 of angle ∆=%.2f°' %angle)
    plt.xlabel('y[um]')
    plt.ylabel('z[um]')

    print("footprint values of optical axis: xoe1=%f, yoe1=%f, zoe2=%f, yoe2=%f" %(xoe1[0]*1e6,yoe1[0]*1e6,zoe2[0]*1e6,yoe2[0]*1e6))

    #plot_scatter(yoe1*1e6, xoe1*1e6, title='footprint on oe1', xtitle='y[um]', ytitle='x[um]')
    #plot_scatter(yoe2*1e6, zoe2*1e6, title='footprint on oe2', xtitle='y[um]', ytitle='z[um]')

def FWHM(beam, Nn):


    Mm = 25
    x = [0] * Mm
    y = [np.ones(Nn)] * Nn

    y2 = [0] * Nn

    for n in range(Nn):

        if beam[n] != 'pollo':

            x = [0] * Mm
            y[n] = [0] * Mm

            print('iteration %d' % n)

            for i in range(Mm):

                x1 = -80 * 1e-6 + (80 + 80) * 1e-6 * i / Mm
                x2 = -80 * 1e-6 + (80 + 80) * 1e-6 * (i + 1) / Mm
                x[i] = -80 * 1e-6 + (80 + 80) * 1e-6 * i / Mm

                for ii in range(beam[n].N):
                    if beam[n].vx[ii] < x2 and beam[n].vx[ii] > x1:
                        y[n][i] = y[n][i] + 1
                if n <4:
                    y[n][i] = y[n][i]*1e-5*10   #+ 2 * n
                else:
                    y[n][i] = y[n][i]*1e-5*10   #+ 2 * (n+1)

    fwhm = np.ones(Nn)
    alphas = np.ones(Nn)
    i = 0
    count_zero = 0
    ciao = 1

    for cc in range(1,Nn):
        if beam[cc] != 'pollo':
            maks = max(beam[cc].vx)
            mins = min(beam[cc].vx)
            d = maks - mins

            print("d=%f" % d)

            fwhm[i] = 2 * np.sqrt(2 * np.log(2)) * d / 10
            alphas[i] = - round((alpha_initial + alpha_step  * cc) * 1e3, 5)

            print("iteration %d" %cc)
            print(alphas[i])

            if alphas[i] < 0 and ciao !=0:
                fwhm[i+1] = 2 * np.sqrt(2 * np.log(2)) * d / 10
                alphas[i+1] = - round((alpha_initial + alpha_step * cc) * 1e3, 5)

                max(beam[0].vx) - min(beam[0].vx)
                fwhm[i] = 2 * np.sqrt(2 * np.log(2)) * d / 10
                alphas[i] = 0.0
                i = i + 1
                ciao = 0

            i = i + 1
        else:
            count_zero=1




    #fwhm[0] = fwhm[1]
    #alphas[0] = alphas[1]

    #alphas[0] = 0.0
    #print("prova")
    #d = max(beam[0].vx) - min(beam[0].vx)
    #fwhm[0] = 2 * np.sqrt(2 * np.log(2)) * d / 7
    #print(fwhm[0])

    alphas[Nn-1] = alphas[Nn-2]
    fwhm[Nn-1] = fwhm[Nn-2]

    print(fwhm)
    print(alphas)

    plt.figure()
    plt.plot(alphas, fwhm*1e6, marker='o', markeredgecolor='r')
    plt.xlabel("∆ in millidegree")
    plt.ylabel("fwhm [urad]")
    plt.title(title)



    if save is True:

        f = h5py.File('dati/FWHM/' + title + 'h5', 'w')
        f['alpha'] = alphas
        f['FWHM'] = fwhm
        f.close()




    d = round(Nn/6)
    salpha = [0] * Nn
    ssalpha = [0] * Nn


    plt.figure()
    plt.plot(x, y[0])
    plt.plot(x, y[1])
    plt.plot(x, y[1*d])
    plt.plot(x, y[2*d])
    plt.plot(x, y[3*d])
    plt.plot(x, y[4*d])
    plt.plot(x, y[5*d])
    plt.legend(['0.0', str(alpha_initial), str(alpha_initial + alpha_step *1 *d), str(alpha_initial + alpha_step * 2*d), str(alpha_initial + alpha_step * 3*d), str(alpha_initial + alpha_step * 4*d), str(alpha_initial + alpha_step * 5*d)])


def FWHM2(beam):

    n, bins, patches = plt.hist(beam.vx, 900, normed=1, histtype='step')
    plt.close()

    (mu, sigma) = norm.fit(beam.vx)
    y = mlab.normpdf(bins, mu, sigma)

    return y, bins, sigma


def plot_FWHM2(beam, Nn):

    y = [0] * Nn
    x = [0] * Nn
    sigma = np.ones(Nn)
    alpha = np.ones(Nn)
    fwhm = np.ones(Nn)

    for i in range (Nn-1):

        print("Iteration %d" %i)

        if beam[i+1] != 'pollo':
            y[i], x[i], sigma[i] = FWHM2(beam[i+1])
            sigma[i] *= 1e6
            fwhm[i] = 2 * np.sqrt(2 * np.log(2)) * sigma[i]
            alphas = - round((alpha_initial + alpha_step  * (i)) * 1e3, 5)
            if alphas <= 0:
                alpha[i] = - round((alpha_initial + alpha_step  * (i+1)) * 1e3, 5)
            else:
                alpha[i] = - round((alpha_initial + alpha_step  * (i)) * 1e3, 5)
        else:
            y[i], x[i], sigma[i] = FWHM2(beam[0])
            sigma[i] *= 1e6
            fwhm[i] = 2 * np.sqrt(2 * np.log(2)) * sigma[i]
            alpha[i] = 0.0




    alpha[Nn-1] = alpha[Nn-2]
    fwhm[Nn-1] = fwhm[Nn-2]

    print(alpha)
    print(fwhm)

    i1 = np.where(fwhm == min(fwhm))

    print(i1[0])

    print("minimum angle %f" %alpha[i1[0][0]])
    print("minimum fwhm %f" %fwhm[i1[0][0]])

    plt.figure()
    plt.plot(alpha, fwhm, marker='o', markeredgecolor='r')

def plot_histogram(beam,Nn):

    d = round(Nn/6)
    salpha = [0.0] * Nn


    salpha[0] = 0.0
    salpha[1] = -alpha_initial

    for i in range (1,6):
        salpha[i+1] = - round(alpha_initial + alpha_step * i * d, 7)


    if x_shift == True:

        beam[1].vx += 10e-6
        for i in range (1,6):
            beam[d*i].vx += 10e-6 * (i+1)



    plt.figure()
    plt.hist(beam[0].vx * 1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(beam[1].vx * 1e6, bins=900, normed=1, histtype='step', color='orange')
    plt.hist(beam[1*d].vx * 1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(beam[2*d].vx * 1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(beam[3*d].vx * 1e6, bins=900, normed=1, histtype='step', color='darkviolet')
    plt.hist(beam[4*d].vx * 1e6, bins=900, normed=1, histtype='step', color='saddlebrown')
    plt.hist(beam[5*d].vx * 1e6, bins=900, normed=1, histtype='step', color='lightpink')
    #plt.hist(beam[6*d].vx * 1e6, bins=900, normed=1, histtype='step', color='pink')
    plt.legend(salpha)
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")

    if show == True:
        plt.show()


def small_plot_histogram(beam, Nn):

    salpha = [0] * Nn

    for i in range(Nn):
        salpha[i] = round(alpha_initial + alpha_step * i, 3)

    plt.figure()
    plt.hist(beam[0].vx * 1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(beam[1].vx * 1e6, bins=900, normed=1, histtype='step', color='orange')
    plt.hist(beam[2].vx * 1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(beam[3].vx * 1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(beam[4].vx * 1e6, bins=900, normed=1, histtype='step', color='darkviolet')
    plt.hist(beam[5].vx * 1e6, bins=900, normed=1, histtype='step', color='saddlebrown')
    plt.hist(beam[6].vx * 1e6, bins=900, normed=1, histtype='step', color='lightpink')
    # plt.hist(beam[6*d].vx * 1e6, bins=900, normed=1, histtype='step', color='pink')
    plt.legend(salpha)
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")


def big_extract_beam_from_file(Nn):

    beam = [0] * Nn
    alpha = [0] * Nn
    salpha = [0] * Nn
    footprint = [0] * Nn

    beam[0] = import_beam(str(0.0))
    footprint[0] = import_footprint(str(0.0))
    salpha[0] = str(0.0) + ' degree'

    for i in range(1, Nn):
        print("Extract iteration %d" %i)
        alpha[i] = round(alpha_initial + alpha_step * i, 5)
        if abs(alpha[i]) < 1e-13:
            beam[i] = "pollo"
            footprint[i] = "pollo"
        else:
            beam[i] = import_beam(str(alpha[i]))
            footprint[i] = import_footprint(str(alpha[i]))



    return beam, footprint

def small_extract_beam_from_file(Nn):

    beam = [0] * Nn
    alpha = [0] * Nn
    salpha = [0] * Nn
    footprint = [0] * Nn

    for i in range(Nn):
        alpha[i] = round(alpha_initial + alpha_step * i, 3)
        beam[i] = import_beam(str(alpha[i]))
        footprint[i] = import_footprint(str(alpha[i]))


    return beam, footprint


if main == '__Moretti_1__':

    filename = 'M.Moretti  29-08-18 at 15:08:53'
    filename = 'dati/Moretti/' + filename + '.h5'

    f = h5py.File(filename, 'r')

    Nn = 7
    beam = [0]*Nn
    alpha = [0]*Nn
    salpha = [0]*Nn
    footprint = [0]*Nn

    for i in range (Nn):

        alpha[i] = round(-0.03 + 0.01 * i, 3)
        beam[i] = import_beam(str(alpha[i]))
        footprint[i] = import_footprint(str(alpha[i]))

        salpha[i] = str(-alpha[i]) + '    degree'

        beam[i].plot_xpzp()
        plt.title("x'/v' plot with an angle of %s" %salpha[i])
        plt.close()

    plt.figure()
    plt.hist(beam[0].vx*1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(beam[1].vx*1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(beam[2].vx*1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(beam[3].vx*1e6, bins=900, normed=1, histtype='step', color='k')
    plt.hist(beam[4].vx*1e6, bins=900, normed=1, histtype='step', color='m')
    plt.hist(beam[5].vx*1e6, bins=900, normed=1, histtype='step', color='c')
    plt.hist(beam[6].vx*1e6, bins=900, normed=1, histtype='step', color='pink')
    plt.legend(salpha)
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")
    plt.title("1um*1um source with gaussian divergence of FWHM 25urad")

    #for i in range (7):
    #    plot_footprint(footprint[i], alpha[i])


    plot_footprint(footprint[3], 0.0)



    f.close()

    plt.show()

if main == '__Moretti_2__':

    filename = 'M.Moretti  31-08-18 at 17:47:45'
    filename = 'dati/Moretti/' + filename + '.h5'

    f = h5py.File(filename, 'r')

    Nn = 100
    beam = [0]*Nn
    alpha = np.ones(Nn)
    salpha = [0]*Nn
    footprint = [0]*Nn

    for i in range (Nn):

        alpha[i] = round(-0.03 + 0.06 * i / Nn, 7)
        beam[i] = import_beam(str(alpha[i]))
        footprint[i] = import_footprint(str(alpha[i]))

        salpha[i] = str(-alpha[i]) + '    degree'

        beam[i].plot_xpzp()
        plt.title("x'/v' plot with an angle of %s" %salpha[i])
        plt.close()

    print(len(salpha))

    Mm = 1000
    x = [0] * Mm
    y = [0] * Mm

    y2 = [0] * Nn

    for n in range (Nn):

        x = [0] * Mm
        y = [0] * Mm

        print('iteration %d' % n)

        for i in range (Mm):

            x1 = -80*1e-6 + 2 * 80 * 1e-6 * i / Mm
            x2 = -80*1e-6 + 2 * 80 * 1e-6 * ( i + 1 ) / Mm
            x[i] =  -80 + 2 * 80 * i / Mm

            for ii in range (beam[n].N):
                if beam[n].vx[ii] < x2 and beam[n].vx[ii] > x1:
                    y[i] = y[i] + 1

            y2[n] = max(y)

    plt.figure()
    plt.plot(alpha, y2)
    plt.title("2*FWHM")

    print(alpha)

    plt.figure()
    indice1 = np.where(alpha == -0.0048)
    c1= indice1[0][0]
    print(c1)
    plt.hist(beam[c1].vx)
    plt.title(salpha[c1])




    plt.figure()
    indice2 = np.where(alpha == 0.0162)
    c2= indice2[0][0]
    print(c2)
    plt.hist(beam[c2].vx)
    plt.title(salpha[c2])

    ss = c2-c1
    fwhm = np.ones(ss)
    alphas = np.ones(ss)
    i=0


    for cc in range (c1, c2):


        maks = max(beam[cc].vx)
        mins = min(beam[cc].vx)
        d = maks - mins


        print("d=%f" %d)

        fwhm[i] = 2*np.sqrt(2*np.log(2))*d/8
        alphas[i] = -alpha[cc]*1e3
        i=i+1

    print(fwhm)
    print(alphas)

    plt.figure()
    plt.plot(alphas, fwhm*1e6)
    plt.xlabel("∆ in millidegree")
    plt.ylabel("fwhm [urad]")
    plt.title("FWHM")

    plt.show()

    indic = np.where(fwhm == np.min(fwhm))

    print(indic)
    print("The perfect deviation is at %f millidegree" %alphas[indic])

    plt.show()

if main == '__Moretti_3__':


    filename = 'M.Moretti  30-08-18 at 15:54:37'
    filename = 'dati/Moretti/' + filename + '.h5'

    f = h5py.File(filename, 'r')

    Nn = 7
    beam = [0]*Nn
    alpha = [0]*Nn
    salpha = [0]*Nn
    footprint = [0]*Nn

    for i in range (Nn):

        alpha[i] = round(-0.03 + 0.01 * i, 3)
        beam[i] = import_beam(str(alpha[i]))
        footprint[i] = import_footprint(str(alpha[i]))

        salpha[i] = str(-alpha[i]) + '    degree'

        beam[i].plot_xpzp()
        plt.title("x'/v' plot with an angle of %s" %salpha[i])
        plt.close()

    alphai = 0.0078
    beami = import_beam(str(alphai))
    footprinti = import_footprint(str(alphai))

    salphai = str(alphai) + '    degree'

    beami.plot_xpzp()
    plt.title("x'/v' plot with an angle of %s" % salphai)
    plt.close()



    plt.figure()
    plt.hist(beam[0].vx*1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(beam[1].vx*1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(beam[2].vx*1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(beam[3].vx*1e6, bins=900, normed=1, histtype='step', color='k')
    plt.hist(beami.vx*1e6, bins=900, normed=1, histtype='step',   color='y')
    plt.hist(beam[4].vx*1e6, bins=900, normed=1, histtype='step', color='m')
    plt.hist(beam[5].vx*1e6, bins=900, normed=1, histtype='step', color='c')
    plt.hist(beam[6].vx*1e6, bins=900, normed=1, histtype='step', color='pink')
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")
    plt.title("1um*1um source with gaussian divergence of FWHM 25urad")
    plt.legend(['0.03 °', '0.02°', '0.01°', '0.0°', '0.0078°', '-0.01°', '-0.02°', '-0.03°'])

    #for i in range (7):
    #    plot_footprint(footprint[i], alpha[i])


    plot_footprint(footprint[3], 0.0)



    f.close()

    plt.show()


if main == '__Moretti_4__':


    filename = 'M.Moretti  31-08-18 at 18:32:35'
    filename = 'dati/Moretti/' + filename + '.h5'

    f = h5py.File(filename, 'r')

    Nn = 7
    beam = [0]*Nn
    alpha = [0]*Nn
    salpha = [0]*Nn
    footprint = [0]*Nn

    for i in range (Nn):

        alpha[i] = round(-0.03 + 0.01 * i, 3)
        beam[i] = import_beam(str(alpha[i]))
        footprint[i] = import_footprint(str(alpha[i]))

        salpha[i] = str(-alpha[i]) + '    degree'

        beam[i].plot_xpzp()
        plt.title("x'/v' plot with an angle of %s" %salpha[i])
        plt.close()

    alphai = 9.600000
    beami = import_beam(str(alphai))
    footprinti = import_footprint(str(alphai))

    salphai = str(alphai) + '    degree'

    beami.plot_xpzp()
    plt.title("x'/v' plot with an angle of %s" % salphai)
    plt.close()

    #for i in range (7):
    #    plot_footprint(footprint[i], alpha[i])


    #plot_footprint(footprint[3], 0.0)

    #plt.figure()
    #plt.plot(beam[3].x*1e6, beam[3].z*1e6, color='r', marker='.',markersize=0.1, linestyle='None')
    #plt.xlabel('x[um]')
    #plt.ylabel('z[um]')
    #plt.title('Image size of the beam')


    plt.figure()
    plt.plot(beam[3].vx*1e6, beam[3].vz*1e6, color='b', marker='.',markersize=0.1, linestyle='None')
    plt.xlabel("x'[urad]")
    plt.ylabel("z'[urad]")
    plt.title('Image divergence of the beam')

    plt.show()

    Mm = 10
    x = [0] * Mm
    y = [np.ones(Nn)] * Mm

    y2 = [0] * Nn

    for n in range(Nn):

        x = [0] * Mm
        y[n] = [0] * Mm

        print('iteration %d' % n)

        for i in range(Mm):

            x1 = -80 * 1e-6 + 2 * 80 * 1e-6 * i / Mm
            x2 = -80 * 1e-6 + 2 * 80 * 1e-6 * (i + 1) / Mm
            x[i] = -80 + 2 * 80 * i / Mm

            for ii in range(beam[n].N):
                if beam[n].vx[ii] < x2 and beam[n].vx[ii] > x1:
                    y[n][i] = y[n][i] + 1
            if n <4:
                y[n][i] = y[n][i]*1e-5*10 + 2 * n
            else:
                y[n][i] = y[n][i]*1e-5*10 + 2 * (n+1)


        f.close()





    plt.figure()
    plt.hist(beam[0].vx*1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(beam[1].vx*1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(beam[2].vx*1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(beam[3].vx*1e6, bins=900, normed=1, histtype='step', color='k')
    plt.hist(beami.vx*1e6, bins=900, normed=1, histtype='step',   color='y')
    plt.hist(beam[4].vx*1e6, bins=900, normed=1, histtype='step', color='m')
    plt.hist(beam[5].vx*1e6, bins=900, normed=1, histtype='step', color='c')
    plt.hist(beam[6].vx*1e6, bins=900, normed=1, histtype='step', color='pink')
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")
    plt.title("1um*1um source with gaussian divergence of FWHM 25urad")
    plt.legend(['∆=0.03 °', '∆=0.02°', '∆=0.01°', '∆=0.0°', '∆=-0.0096°', '∆=-0.01°', '∆=-0.02°', '∆=-0.03°'])



    ym = np.ones(Mm)
    ym *= 0.

    for i in range(Mm):

        x1 = -80 * 1e-6 + 2 * 80 * 1e-6 * i / Mm
        x2 = -80 * 1e-6 + 2 * 80 * 1e-6 * (i + 1) / Mm
        x[i] = -80 + 2 * 80 * i / Mm

        for ii in range(beami.N):
            if beami.vx[ii] < x2 and beami.vx[ii] > x1:
                ym[i] = ym[i] + 1
        ym[i] = ym[i] * 1e-5 *10 + 2 * 4





    plt.figure(figsize=(5,10))
    plt.plot(x, y[0], color = 'r')
    plt.plot(x, y[1], color = 'b')
    plt.plot(x, y[2], color = 'g')
    plt.plot(x, y[3], color = 'k')
    plt.plot(x, ym  , color = 'pink')
    plt.plot(x, y[4], color = 'y')
    plt.plot(x, y[5], color = 'm')
    plt.plot(x, y[6], color = 'c')
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")
    plt.title("1um*1um source with gaussian divergence of FWHM 25urad")
    plt.legend(['∆=0.03 °', '∆=0.02°', '∆=0.01°', '∆=0.0°', '∆=-0.0096°', '∆=-0.01°', '∆=-0.02°', '∆=-0.03°'])


    plt.figure(figsize=(5,10))
    plt.plot(x, ym)


    plt.show()



if main == '__Moretti_5__':

#    filename = 'M.Moretti  05-09-18 at 09:54:42'
#    filename = 'dati/Moretti/' + filename + '.h5'


    f = h5py.File(filename, 'r')

    if dim == 'big':
        Nn = 80
        beam, footprint = big_extract_beam_from_file(Nn=Nn)
        plot_histogram(beam,Nn)
        f.close()
        plot_footprint(footprint[0], 0.0)
        FWHM(beam,Nn)

        beam[0].plot_xpzp()

    elif dim == 'small':
        Nn = 7
        beam, footprint = small_extract_beam_from_file(Nn)
        plot_footprint(footprint[3], 0.0)
        small_plot_histogram(beam, Nn)

        beam[3].plot_xpzp()

        f.close()




    plt.show()



if main == '__Moretti_7__':

    f = h5py.File(filename, 'r')


    Nn = 80
    beam, footprint = big_extract_beam_from_file(Nn)
    plot_histogram(beam, Nn)
    f.close()

    plot_FWHM2(beam, Nn)
    plot_footprint(footprint[0], 0.0)

    plot_histogram(beam, Nn)

    beam[0].plot_xz()
    plt.title("beam at the end of the Montel")
    beam[0].plot_xpzp()
    plt.title("divergence at the end of the Montel")


    plt.show()