import healpy as hp
import numpy as np

from scipy.optimize import curve_fit
from scipy import stats

def maskMap(m, decmin, decmax):

    degree = np.pi / 180.
    npix  = len(m)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    thetaMin = (90 - decmin) * degree
    thetaMax = (90 - decmax) * degree
    thetaCut = (theta <= thetaMin) * (theta >= thetaMax)

    new_map = np.copy(m)
    new_map[np.logical_not(thetaCut)] = hp.UNSEEN

    return new_map

def multipole2dfunc(l,nside):

    def multipole2d(ipix, *p):
        k = np.pi/180         # Wavenumber (assumes x in degree)

        ra,dec = hp.pix2ang(nside, ipix,lonlat=True)

        horizontal = sum([p[2*i] * np.power(np.cos(k*dec),i+1) * np.cos((i+1)*k*ra - p[2*i+1]) for i in range(l)])
        return horizontal

    return multipole2d

def multipoleFit2d(relint, relerr, l, nside_out=None, decmin=-90., decmax=-35., no_err=False):

    # Mask maps
    relint = maskMap(relint, decmin, decmax)
    relerr = maskMap(relerr, decmin, decmax)

    # I feel like the masked pixels need to be set to 0 with an infinite uncertainty
    # Check that, doesn't seem to matter
    #relint[relint == hp.UNSEEN] = 0
    #relerr[relerr == hp.UNSEEN] = np.inf
    
    variance = np.power(relerr,2)
    if nside_out:
        relint = hp.ud_grade(relint, nside_out)
        variance= hp.ud_grade(variance, nside_out, power=2)
        relerr=np.sqrt(variance)

    npix = relint.size
    nside = hp.npix2nside(npix)
    minpix = hp.ang2pix(nside,0,-30,lonlat=True)
    pixels = range(minpix,npix)
    maxphase=2*np.pi

    # Guess at best fit parameters
    amplitude = 1e-3
    phase     = 1
    p0 = [amplitude, phase] * l

    # Define bounds
    b0 = [0,0] * l
    b1 = [np.inf,maxphase] * l
    
    fitfunc = multipole2dfunc(l,nside)
    if no_err:
        popt, pcov = curve_fit(fitfunc,
            pixels,
            relint[pixels],
            p0 = p0
            #bounds=(b0,b1), #for some reason the combination of b0 and p0 doen't work
            #sigma=relerr[pixels], absolute_sigma=True
            )
    else:
        popt, pcov = curve_fit(fitfunc,
            pixels,
            relint[pixels],
            p0,
            #bounds=(b0,b1), #for some reason the combination of b0 and p0 doen't work
            sigma=relerr[pixels], absolute_sigma=True
            )
    chi2 = sum((relint[pixels] - fitfunc(pixels, *popt))**2 / relerr[pixels]**2)
    ndof = len(pixels) - popt.size
    pvalue = 1 - stats.chi2.cdf(chi2, ndof)
    perr = np.sqrt(np.diag(pcov))

    return popt, perr, chi2, ndof, pvalue


# Unweighted average consistent with 6-year publication
def returnRI(relint_map, relerr_map, nbins=24, decmin=-90., decmax=-35.):

    # Mask maps
    relint_map = maskMap(relint_map, decmin, decmax)
    relerr_map = maskMap(relerr_map, decmin, decmax)

    # Setup right-ascension bins
    degree = np.pi / 180
    ramin, ramax = 0, 360 * degree
    rabins = np.linspace(ramin, ramax, nbins+1)

    # Calculate phi for each pixel
    npix  = len(relint_map)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    # Bin in right ascension
    phiBins = np.digitize(phi, rabins) - 1
    # UNSEEN cut
    cut = (relint_map != hp.UNSEEN)

    ri, sigmay = np.zeros((2,nbins))
    for i in range(nbins):
        phiCut = (phiBins == i)
        c0 = cut * phiCut
        ri[i] = np.mean(relint_map[c0])
        # Error result from propagation of uncertainty on unweighted average
        sigmay[i] = np.sqrt(np.sum(relerr_map[c0]**2))/c0.sum()

    dx = (ramax - ramin)/(2*nbins)
    ra = np.linspace(ramin+dx, ramax-dx, nbins) / degree
    sigmax = dx * np.ones(nbins) / degree

    return (ra, ri, sigmax, sigmay)


# Cosine function with fixed base wavelength of 360 degrees
def multipole(x, *p):
    k = 2*np.pi/360         # Wavenumber (assumes x in degrees)
    l = int(len(p) / 2)    # Multipole number of fit
    return sum([p[2*i] * np.cos((i+1)*k*x - p[2*i+1]) for i in range(l)])
    

# Best-fit parameters for multipole
def multipoleFit1D(x, y, l, sigmay):

    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    p0 = [amplitude, phase] * l

    # Do best fit
    popt, pcov = curve_fit(multipole, x, y, p0,
                           sigma=sigmay, absolute_sigma=True)
    chi2 = sum((y - multipole(x, *popt))**2 / sigmay**2)

    return popt, pcov, chi2
