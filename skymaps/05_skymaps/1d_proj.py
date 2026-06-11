#!/usr/bin/env python

from __future__ import print_function


import numpy as np
import healpy as hp
import os, optparse, re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.container import ErrorbarContainer
from scipy.optimize import curve_fit
from scipy import stats
from scipy.integrate import quad
from mapFunctions import getMap

#rc('xtick', labelsize=18)
#rc('ytick', labelsize=18)

#-------------------------------------------------------

def getRelErr(f):
    return np.sqrt(hp.read_map(f, 0, verbose=False)**2)

def getRelInt(f, rel_err=None, **opts):
    
    relint = getMap([f], mapName='relint', **opts)
    if rel_err:
        relerr = getRelErr(rel_err)
    else:
        relerr = getMap([f], mapName='relerr', **opts)

    # Setup right-ascension bins
    degree = np.pi / 180
    ramin = opts['ramin'] * degree
    ramax = opts['ramax'] * degree
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)

    # Calculate phi for each pixel
    npix  = len(relint)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    phiBins = np.digitize(phi, rabins) - 1


    # Treat masked pixels as zeroes for weighting purposes
    cut = (relint != hp.UNSEEN)

    ri, sigmay = np.zeros((2, opts['nbins']))

    for i in range(opts['nbins']):
        phiCut = (phiBins == i)
        c0 = cut * phiCut
        ri[i] = np.mean(relint[c0])
        # Error result from propagation of uncertainty on unweighted average
        sigmay[i] = np.sqrt(np.sum(relerr[c0]**2))/c0.sum()

    dx = (ramax - ramin)/(2*opts['nbins'])
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / degree
    sigmax = dx * np.ones(opts['nbins']) / degree

    return (ra, ri, sigmax, sigmay)

# Flat line
def flatLine(x, *p):
    return 0*x + p[0]

# Best-fit parameters for flat line
def lineFit(x, y, sigmay):
    popt, pcov = curve_fit(flatLine, x, y, [0], 
                           sigma=sigmay, absolute_sigma=True)
    chi2 = sum((y - flatLine(x, *popt))**2 / sigmay**2)
    return popt, pcov, chi2


# Cosine function with fixed base wavelength of 360 degrees
def multipole(x, *p):
    k = 2*np.pi/360         # Wavenumber (assumes x in degrees)
    l = int(len(p) / 2)    # Multipole number of fit
    return sum([p[2*i] * np.cos((i+1)*k*x - p[2*i+1]) for i in range(l)])

# Best-fit parameters for multipole
def multipoleFit(x, y, l, sigmay):

    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    p0 = [amplitude, phase] * l

    # Do best fit
    popt, pcov = curve_fit(multipole, x, y, p0, 
                           sigma=sigmay, absolute_sigma=True)
    chi2 = sum((y - multipole(x, *popt))**2 / sigmay**2)

    return popt, pcov, chi2

# Gaussian function
def gaussian(x, *p):
    numerator = x-p[1]
    denominator = np.sqrt(2)*p[2]
    return p[0] * np.exp(-(numerator/denominator)**2) + p[3]

def gaussian_for_plotting(x, *p):
    x0 = np.remainder((x + 90), 360)
    numerator = x0 - (p[1] + 90)
    denominator = np.sqrt(2)*p[2]
    return p[0] * np.exp(-(numerator/denominator)**2) + p[3]

# Function to fit the Gaussian and adjust the baseline
def fit_gaussian_with_zero_integral(x_data, y_data, x_min, x_max, sigmay):
    # Initial guess for amplitude, mean, stddev, and baseline
    amplitude = -3.67
    ra_cod = 94.1
    width = 65.3
    b = 1.58
    p0 = [amplitude, ra_cod, width, b]
    #x0 = np.remainder((x_data+ 90), 360)
    # Wrapper for curve_fit that enforces the integral constraint
    def constrained_gaussian(x, amplitude, ra_cod, width):

        # Find the baseline that makes the integral zero
        baseline = -quad(lambda x: amplitude * np.exp(-((x - ra_cod) / (np.sqrt(2)*width))**2), (x_min-90), (x_max-90))[0] / (x_max - x_min)
   
        p = [amplitude, ra_cod, width, baseline]
        return gaussian(x, *p)
    
    # Fit the data
    popt, pcov = curve_fit(constrained_gaussian, x_data, y_data, p0=p0[:-1], sigma = sigmay, absolute_sigma = True)
    # Calculate the final baseline
    baseline = -quad(lambda x: popt[0] * np.exp(-((x - popt[1]) / (np.sqrt(2) * popt[2]))**2), (x_min-90), (x_max-90))[0] / (x_max - x_min)
    popt = np.append(popt, baseline)
   # popt[1] = np.remainder((popt[1] - 90), 360) 
    
    # Adding the baseline term to the covariance matrix
    amplitude, mean, stddev = popt[:-1]
    
    def integral_derivative_wrt_amplitude(x):
        return np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))
    
    def integral_derivative_wrt_mean(x):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2)) * (x - mean) / (stddev ** 2)
    
    def integral_derivative_wrt_stddev(x):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2)) * ((x - mean) ** 2) / (stddev ** 3)
    
    deriv_integral_amplitude = quad(integral_derivative_wrt_amplitude, x_min, x_max)[0]
    deriv_integral_mean = quad(integral_derivative_wrt_mean, x_min, x_max)[0]
    deriv_integral_stddev = quad(integral_derivative_wrt_stddev, x_min, x_max)[0]
    
    variance_baseline = (
        (deriv_integral_amplitude / (x_max - x_min))**2 * pcov[0, 0] +
        (deriv_integral_mean / (x_max - x_min))**2 * pcov[1, 1] +
        (deriv_integral_stddev / (x_max - x_min))**2 * pcov[2, 2] +
        2 * (deriv_integral_amplitude * deriv_integral_mean) / (x_max - x_min)**2 * pcov[0, 1] +
        2 * (deriv_integral_amplitude * deriv_integral_stddev) / (x_max - x_min)**2 * pcov[0, 2] +
        2 * (deriv_integral_mean * deriv_integral_stddev) / (x_max - x_min)**2 * pcov[1, 2]
    )
    
    # Expand the covariance matrix to include the baseline term
    pcov_expanded = np.zeros((4, 4))
    pcov_expanded[:-1, :-1] = pcov
    pcov_expanded[-1, -1] = variance_baseline
    
    
    # Calculate the chi-squared value
    chi_squared = np.sum((y_data - gaussian_for_plotting(x_data, *popt)) ** 2 / sigmay**2)
    
    return popt, pcov_expanded, chi_squared
    
# Best-fit parameters for gaussian
def simple_gaussian_fit(x, y, sigmay):

    # Guess at best fit parameters
    amplitude = -3
    ra_cod = 88
    width = 43.1
    b = 1
    p0 = [amplitude, ra_cod, width, b]

    # Do best fit
    popt, pcov = curve_fit(gaussian, x, y, p0, sigma = sigmay, absolute_sigma = True)
    chi2 = sum((y - gaussian(x, *popt))**2/sigmay**2)

    return popt, pcov, chi2

    #03282023 ADD SYSTEMATIC ERROR ESTIMATION
def systematicErr(err_ri):
    
    # Load data
    # High Energy Anti-Sidereal
    #as_ra, as_ri, as_sigmay = np.loadtxt('./recofits/10YAntiSiderealHigh.txt', delimiter=',', unpack=True)

    #Make histogram?
    #count, bin_edge = np.histogram(y)
    #bin_centre = (bin_edge[:-1] + bin_edge[1:])/2
    #error = np.sqrt(count)

    #plt.errorbar(bin_centre, count, xerr=0, yerr=error, marker='.')
    #plt.savefig('./histogramas.png')

    #print mean and st. dev. on histogram
    mean = np.mean(err_ri)
    stdev = np.std(err_ri)
    #stdev = np.sqrt(np.mean(err_ri**2))
    #print(mean)
    #print(stdev)
    return stdev

    
##=======================================================================##

if __name__ == "__main__":

    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = optparse.OptionParser(usage)
    parser.add_option("-r", "--ramin", dest="ramin", type=float,
            default=0, help="minimum RA")
    parser.add_option("-R", "--ramax", dest="ramax", type=float,
            default=360, help="maximum RA")
    parser.add_option("--rimin", dest="rimin", type=float,
            default=.001, help="minimum relative intensity")
    parser.add_option("--rimax", dest="rimax", type=float,
            default=.001, help="maximum relative intensity")
    parser.add_option("-D", "--decmax", dest="decmax", type=float,
            help="maximum Dec")
    parser.add_option("-d", "--decmin", dest="decmin", type=float,
            help="minimum Dec")
    parser.add_option("-n", "--nbins", dest="nbins", type=int,
            default=24, help="number of bins")
    parser.add_option("-z","--zero", action="store_true", dest="zeroline",
            default=False, help="Draw zero line")
    parser.add_option("-f","--flipra", action="store_true", dest="flipra",
            default=False, help="Flips RA in x axis")
    parser.add_option("-o", "--output", dest="output", default=None,
            help="Output image file name")
    parser.add_option("--multi", dest='multi', type=int,
            default=None, help='Use multipole subtraction')
    parser.add_option('--multiErr', dest='multiErr',
            default=False, action='store_true',
            help='Use amplitude of dipole fit to calculate sys error')
    parser.add_option('--fit', dest='fit',
            type=int,
            help='Show best-fit multipole on plot (1=dipole, 2=quad, etc.)')
    parser.add_option('--offset', dest='offset',
            default=False, action='store_true',
            help='Offset points to avoid overlap')
    parser.add_option('--split', dest='split',
            default=False, action='store_true',
            help='Split input, sharing x-axis')
    parser.add_option("--labels", dest='labels',
            help='Custom label options built-in [configs, method]')
    parser.add_option("-v","--verbose", action="store_true", dest='verbose',
            default=False, help='Optional additional output')
    parser.add_option('--full', dest='full',
            default=False, action='store_true',
            help='Show average behavior of full time range')
    parser.add_option('--gaussian', dest='gaussian',
            default=False, action='store_true',
            help='Fit with a Gaussian function')
    parser.add_option('--title', dest='title',    
            default = None,
            help='add plot title')
    parser.add_option('--sci', dest='sci_notation',
            default=False, action='store_true',
            help='Show probability of fit in scientific notation')
    parser.add_option('--syserr', dest='sys_err',
            default=None,
            help='Produce histogram of antisidereal data points to estimate systematic error')
    parser.add_option('--totalerr', dest='total_err', default=False, action='store_true', help='Show grey as total vs sys err')
    parser.add_option('--relerr', dest="rel_err",
            default=None,
            help='File with relative errors')
    
    # NOTE: I'd love to clean this up, but I think I just need to add to it. 
    # Need: scale (multiply by RI and show up on axis)
    parser.add_option('-S', '--scale', dest='scale',
            type=int, default=0,
            help='Exponential scale for multiplying y-axis')
    parser.add_option('-L', '--legend', dest='legend',
            default=False, action='store_true',
            help='Display plot legend')
    parser.add_option('--flat', dest='flat',
            default=False, action='store_true',
            help='Show fit parameters for a flat line')
    parser.add_option('-m', '--min', dest='min',
            type=float,
            help='Set plot minimum value')
    parser.add_option('-M', '--max', dest='max',
            type=float,
            help='Set plot maximum value')

    options, args = parser.parse_args()
    opts = vars(options).copy()

    # Default masking behavior
    if not options.decmax:
        opts['mask'] = True

    if options.verbose:
        for key in sorted(opts.keys()):
            print(' --%s %s' % (key, opts[key]))

    rel_err = options.rel_err
    
    # Setup dictionaries for default values
    methods = {'sid':'sidereal', 'solar':'solar', 'anti':'anti-sidereal'}
    methods['ext'] = 'extended-sidereal'
    errDict = {'sid':'anti', 'solar':'ext'}
    p = {}      # Storage for best-fit parameters

    # Plotting setup
    figdims = (8,6)
    if opts['offset']:
        figdims = (17,6)
    axs = None
    if opts['split']:
        fig, axs = plt.subplots(1 + additional_graphs, sharex=True, sharey=True, figsize=figdims)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    else:
        fig = plt.figure(figsize=figdims)
        ax = fig.add_subplot(111)

    # Setup for referee adjustments
    ramin = opts['ramin']
    ramax = opts['ramax']
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)

    # Systematic error bar setup
    eb = {'edgecolor':None, 'linewidth':0, 'linestyle':None}

    # Y-axis scale
    scale = 10**opts['scale']
    
    for i in range(0, len(args)):
        if (i < len(args)):
            f = args[i]
            # do f stuff - label, err file, ra etc
            # Labeling
            basename = os.path.basename(f)[:-5]
            label = basename
            # Specialty labeling based on configuration or time frame
            if options.labels == 'configs':
                label = re.findall('IC86-\\d{4}', f)[-1]
            # Find time frame, stripping leading/trailing characters
            #method = re.findall(...
            method = "sidereal"
            if options.labels == 'method':
                label = method

            # Possible file for systematic uncertainties
            errfile = ''
            if method in errDict.keys():
                errfile = f.replace(method, errDict[method])

            # Get RA, relint, and errors
            ra, ri, sigmax, sigmay = getRelInt(f, **opts)
            ri *= scale
            sigmay *= scale
            

            # Offset
            if opts['offset']:
                npoints = len(args) + additional_graphs
                dx = (rabins[1] - rabins[0])/float(npoints+1)
                ra = rabins[:-1] + (npoints-i)*dx

        # SYSTEMATIC ERROR MARCH 28 2023
        if options.sys_err != None:
            if options.sys_err.lower() == 'high sid':
                f = './recofits/10YAntiSiderealHigh.txt'
            if options.sys_err.lower() == 'low sid':
                f = './recofits/10YAntiSiderealLow.txt'
            if options.sys_err.lower() == 'high sol':
                f = './recofits/10YExtSiderealHigh.txt'
            if options.sys_err.lower() == 'low sol':
                f = './recofits/10YExtSiderealLow.txt'
            #Plot slightly translucent graph with systematic error bars

            #Set up txt file by running first on antisidereal.
            #with open(f, 'w') as file:
            #    for i in range(len(ra)):
            #        print(str(ra[i]) + "," + str(ri[i]) + "," + str(sigmay[i]), file=file)
            
            err_ra, err_ri, err_sigmay = np.loadtxt(f, delimiter=',', unpack=True)
            sysErr = systematicErr(err_ri)

            #systematicErr or 2*systematicErr and  total Err = stat + sys, or just sys?
            if options.total_err:
                totalErr = sigmay + sysErr
            else:
                totalErr = np.full_like(sigmay, sysErr)

            l = ax.errorbar(ra, ri, xerr=0*sigmax, yerr=totalErr,
                    marker='.', fmt='.',
                    ecolor='darkgrey', linewidth=2, capsize=0, mew=0)
            
        if opts['offset']:
            print(rabins)
            npoints = len(args)
            dx = (rabins[1] - rabins[0])/float(npoints+1)
            ra = rabins[:-1] + (npoints-i)*dx

        # Select plotting axis
        ax = axs[i] if opts['split'] else ax
        
        # Optionally adjust y-axis limits
        if options.min != None and options.max != None:
            ax.set_ylim(options.min, options.max)

        # Find time frame, stripping leading/trailing characters
        #method = re.findall('|'.join([f'_{k}[_|\.]' for k in methods]), f)[-1]
        #print(len(ra), len(ri), len(sigmax), len(sigmay))
        # Plot 1-D relative intensity projection
        l = ax.errorbar(ra, ri, xerr=0*sigmax, yerr=sigmay,
                    marker='.', fmt='.',
                    capsize=4, label=label, linewidth=2, markersize=8, mew=0)
        if options.verbose:
            print("sigma: ", basename, np.std(ri))
            print("max: ", basename, np.max(np.abs(ri)))
        
        # Use Systematic Error bar to find fit?
        if options.sys_err:
            sigmay = totalErr

        # Optional best fit(s)
        if options.fit != None:

            # Calculate best-fit line and uncertainties
            popt, pcov, chi2 = multipoleFit(ra, ri, options.fit, sigmay)
            eopt = np.sqrt(np.diag(pcov))

            # Plot a smoothed version of the applied fit
            tOpt = {'color':'blue'}
            smooth_x = np.linspace(0, 360, 1000)
            ax.plot(smooth_x, multipole(smooth_x, *popt), **tOpt)

            # Calculate degrees of freedom and p-value
            ndof = ra.size - popt.size
            pvalue = 1 - stats.chi2.cdf(chi2, ndof)

            # Text locations and values
            ymax = ax.get_ylim()[-1]
            delta = ymax/10
            x0 = 100  # x-location (in degrees)
            y0 = 0 if method=='sid' else ymax
            # Adjust best-fit amplitude and phase to positive values
            amp = popt[0]
            phase = popt[1] * 180/np.pi
            if amp < 0:
                phase -= 180
                amp *= -1
            while phase < 0:
                phase += 360
            phase_err = eopt[1] * 180/np.pi

            # Write to plot
            tOpt.update({'fontweight':'bold'})
            info = {}
            info[0] = f'A = {amp:0.2f} $\\pm$ {eopt[0]:0.2f}'
            info[1] = f'$\\phi$ = {phase:0.1f} $\\pm$ {phase_err:0.1f}'
            info[2] = f'$\\chi^2$/NDOF = {chi2:0.0f} $/$ {ndof:0.0f}'
            if method != 'sid':
                if options.sci_notation:
                    info[3] = 'Probability = ' + np.format_float_scientific(pvalue, precision=4)
                else:
                    info[3] = f'Probability = {pvalue:0.4f}'
            for i in info.keys():
                ax.text(x0, y0-(i+1)*delta, info[i], **tOpt)

        if options.flat:

            # Calculate and plot flat line
            popt, pcov, chi2 = lineFit(ra, ri, sigmay)
            eopt = np.sqrt(np.diag(pcov))
            tOpt = {'color':'orange'}
            ax.plot(ra, flatLine(ra, *popt), **tOpt)

            # Calculate degrees of freedom and p-value
            ndof = ra.size - popt.size
            pvalue = 1 - stats.chi2.cdf(chi2, ndof)

            # Text locations and values
            ymax = ax.get_ylim()[-1]
            delta = ymax/10
            x0 = 325   # x-location (in degrees)
            y0 = 0 if method=='sid' else ymax

            # Write to plot
            tOpt.update({'fontweight':'bold'})
            info = {}
            info[0] = f'$\\chi^2$/NDOF = {chi2:0.0f} $/$ {ndof:0.0f}'
            if options.sci_notation:
                info[1] = 'Probability = ' + np.format_float_scientific(pvalue, precision=4)
            else:
                info[1] = f'Probability = {pvalue:0.4f}'
            for i in info.keys():
                ax.text(x0, y0-(i+1)*delta, info[i], **tOpt)
        
        if options.gaussian:

            # Calculate and plot gaussian fit
            popt, pcov, chi2 = fit_gaussian_with_zero_integral(ra, ri, 0, 360, sigmay)
            eopt = np.sqrt(np.diag(pcov))

            # Plot a smoothed version of the applied fit
            tOpt = {'color':'green'}
            smooth_x = np.linspace(0, 360, 1000)
            ax.plot(smooth_x, gaussian_for_plotting(smooth_x, *popt), **tOpt)

            # Calculate degrees of freedom and p-value
            ndof = ra.size - popt.size
            pvalue = 1 - stats.chi2.cdf(chi2, ndof)


            # Text locations and values
            ymax = ax.get_ylim()[-1]
            delta = ymax/10
            x0 = 350    # x-location (in degrees)
            y0 = 0 if method=='sid' else ymax
            # Adjust best-fit amplitude and phase to positive values
            #COPIED FROM MULTIPOLE CHECK IF NECESSARY AND MAKE ADJUSTMENTS

            amp = popt[0]
            ra_cod = popt[1]
            width = popt[2]
            b = popt[3]
            
            # Write to plot
            tOpt.update({'fontweight':'bold'})
            info = {}
            info[0] = f'A = {amp:0.2f} $\\pm$ {eopt[0]:0.2f}'
            info[1] = f'RA COD = {ra_cod:0.1f} $\\pm$ {eopt[1]:0.1f}'
            info[2] = f'WIDTH = {width:0.1f} $\\pm$ {eopt[2]:0.1f}'
            info[3] = f'B = {b:0.2f} $\\pm$ {eopt[3]:0.2f}'
            info[4] = f'$\\chi^2$/NDOF = {chi2:0.0f} $/$ {ndof:0.0f}'
            if method != 'sid':
                if options.sci_notation:
                    info[5] = 'Probability = ' + np.format_float_scientific(pvalue, precision=3)
                else:
                    info[5] = f'Probability = {pvalue:0.4f}'
            for i in info.keys():
                ax.text(x0, y0-(i+1)*delta, info[i], **tOpt)

        # Additional plotting options
        if not os.path.isfile(errfile) and errfile!='':
            print(f'Error file {errfile} not found! Cannot calculate syserr')

        if os.path.isfile(errfile):

            ra_err, ri_err, sigx, sigy = getRelInt(errfile, **opts)
            ri_err *= scale
            sigy *= scale

            # Ideal scenario - fit cosine function and take amplitude
            # Not used because anti- and extended-sidereal not sinusoidal
            if options.multiErr:
                popt, pcov, chi2  = multipoleFit(ra_err, ri_err, 1, sigy)
                syserr = popt[0]
            else:
                #syserr = np.abs(ri_err).max()  # Six-year method
                syserr = np.sqrt(np.mean(ri_err**2))

            box = dx if opts['offset'] else 10
            patches = [mpl.patches.Rectangle([ra[j]-box/2, ri[j]-syserr], box,
                    2*syserr, **eb) for j in range(len(ra))]
            cln = PatchCollection(patches, cmap=mpl.cm.jet,
                    alpha=0.5, facecolor=l[0].get_color())
            ax.add_collection(cln)

        # Add dashed lines to separate right ascension bins
        if opts['offset']:
            for raval in rabins[1:-1]:
                ax.axvline(x=raval, color='k', ls=':')
     
        # Axes labels
        ax.set_xlabel(r"Right Ascension $[^{\circ}]$",fontsize=14)
        if opts['split']:
            ylabel = r'$\Delta N/\langle N \rangle$'
            x0, y0 = -0.04, 0.5
            if opts['scale'] != 0:
                ylabel += fr'$\; (\times 10^{{{-opts["scale"]}}})$'
                x0, y0 = 0.05, 0.5
            fig.text(x0, y0, ylabel, fontsize=14,
                    va='center', rotation='vertical')
        else:
            ylabel = r'$\Delta N/\langle N \rangle$'
            if opts['scale'] != 0:
                ylabel += fr'$\;(\times 10^{{{-opts["scale"]}}})$'
            ax.set_ylabel(ylabel,fontsize=14)
        #ax.grid()

    # Optionally include average of whole dataset
    axs = axs if opts['split'] else [ax]
    for ax in axs:
        if options.full:
            totalMap = re.sub('IC86-\\d{4}_','IC86_', f)
            ra, ri, sigmax, sigmay = getRelInt(totalMap, **opts)
            ri *= scale
            sigmay *= scale
            totalErr = re.sub('IC86-\\d{4}_','IC86_', errfile)
            if os.path.isfile(totalErr):
                print(f'Error file {os.path.basename(errfile)} found.')
                ra_err, ri_err, sigx, sigy = getRelInt(totalErr, **opts)
                ri_err *= scale
                sigy *= scale
                if options.multiErr:
                    popt, pcov, chi2 = multipoleFit(ra_err, ri_err, 1, sigy)
                    syserr = popt[0]
                else:
                    #syserr = np.abs(ri_err).max()  # Six-year method
                    syserr = np.sqrt(np.mean(ri_err**2))
                box = dx * len(args)
                patches = [mpl.patches.Rectangle([ra[j]-box/2, ri[j]-syserr],
                        box, 2*syserr, **eb) for j in range(len(ra))]
                cln = PatchCollection(patches, alpha=0.4, facecolor='gray')
                ax.add_collection(cln)

        # Show zero line
        if options.zeroline:
            xzero = np.arange(0, 360, 1)
            yzero = 0 * xzero
            ax.plot(xzero,yzero,linewidth=1.5,linestyle='--',color='black')

        # Adjust x-axis limits
        ax.set_xlim(options.ramax, options.ramin)
        if options.flipra:
            ax.set_xlim(options.ramin, options.ramax)

        # Show legend
        if opts['legend']:
            leg = ax.legend(loc='lower right')
        
        # Add Title
        if opts['title']:
            plt.title(options.title)

    plt.draw()
    plt.savefig(options.output, dpi=300, bbox_inches='tight')
    

