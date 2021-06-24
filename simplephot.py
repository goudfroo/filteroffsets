#!/usr/bin/env python
#
# Syntax: python simplephot.py <input fits file> <nsigma> <outtab>
#
#     Perform astrometry and simple photometry of point-like sources in input
#     image. Calculates clipped mean of the sky annulus as background value. 
#     Eliminates sources with more than 1 saturated pixel in the source aperture
#     and/or "do not use" pixels from the initial source list prior to
#     performing final photometry.
#
#     Optional input parameters:
#     nsigma (#2): Number of sigmas above the background to use as detection 
#                  limit. Defaulted to 100.
#     outtab (#3): Name of (ASCII-format) output photometry table.
#                  Defaulted to "simplephot.out".
#
import sys
import os
import math
import numpy as np
from photutils import DAOStarFinder
import photutils.detection as detection
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry as apphot
from photutils.centroids import centroid_sources, centroid_com
from photutils.morphology import data_properties
from astropy import units
from astropy import stats
from astropy.io import fits

#------------- DEFINITIONS -----------------------------------------
def find_sources(image, fwhm=1.4, nsigma=100.0, roundness=1.0, sharpness=0.1):
    mean, median, standard_deviation = stats.sigma_clipped_stats(image, sigma=5.0)
    newimage = image - median
    print('Median value subtracted during detection: %.3f\n' % (median))
    print('  Image clipped standard deviation value: %.3f\n' % (standard_deviation))
    daofind = DAOStarFinder(fwhm=fwhm, threshold=nsigma*standard_deviation,
                            roundlo= -roundness, roundhi=roundness,
                            sharplo=sharpness)
    source_list = daofind(newimage) 
    print('Number of sources found by DAOfind: ', len(source_list))
    return source_list

def phot_circ(image, dqarr, source_list, outfile, radii, skip=5., w_ann=5., pixscale=0.0656, gain=1.6):
#    def phot_rect(image, source_list, outfile, w=10.0, h=10.0, skip=2., w_ann=5., pa=0.0):
    xarr = source_list['xcentroid']
    yarr = source_list['ycentroid']
    # Refine the center coordinates from DAOfind using the 2D image moments
    newx, newy = centroid_sources(image, xarr, yarr, box_size=7,
                                        centroid_func=centroid_com)
    coords = np.column_stack((newx,newy))
    # Eliminate stars with > 1 saturated pixels to avoid spurious results
    srcaper = CircularAnnulus(coords, r_in=1, r_out=3)
    srcaper_masks = srcaper.to_mask(method='center')
    satflag = np.zeros((len(newx),), dtype=int)
    i = 0
    for mask in srcaper_masks:
        srcaper_dq = mask.multiply(dqarr)
        srcaper_dq_1d = srcaper_dq[mask.data > 0]
        badpix = np.logical_and(srcaper_dq_1d>2, srcaper_dq_1d<7)
        reallybad = np.where(srcaper_dq_1d == 1)
        totbadpix = len(srcaper_dq_1d[badpix]) + len(srcaper_dq_1d[reallybad])
        if ((len(srcaper_dq_1d[badpix]) > 1) or (len(srcaper_dq_1d[reallybad]) > 0)):
            #print('{} bad pixels in source {}'.format(totbadpix,i))
            satflag[i] = 1
        i = i+1
    goodx = newx[np.where(satflag == 0)]
    goody = newy[np.where(satflag == 0)]
    print('Number of sources without saturated or bad pixels: ', len(goodx))
    print(' ')
    coords = np.column_stack((goodx,goody))
    #
    apertures = [CircularAperture(coords, r=r) for r in radii]
    maxrad = max(radii)
    bckaper = CircularAnnulus(coords, r_in=maxrad+skip, r_out=maxrad+skip+w_ann)
#    aperture = RectangularAperture(coords, w=w, h=h, theta=pa)
#    bckaper = RectangularAnnulus(coords, w_in=w+skip, w_out=w+skip+w_ann,
#                                 h_out=h+skip+w_ann, theta=pa)
    bckaper_masks = bckaper.to_mask(method='center')
    bck_mean = []
    bck_stdev = []
    bck_npix = []
    for mask in bckaper_masks:
        bckaper_data = mask.multiply(image)
        bckaper_data_1d = bckaper_data[mask.data > 0]
        mean_sigclip, _, sig_sigclip = stats.sigma_clipped_stats(bckaper_data_1d, sigma=3.0)
        clipdata = stats.sigma_clip(bckaper_data_1d, sigma=3.0, maxiters=10)
        bck_mean.append(mean_sigclip)
        bck_stdev.append(sig_sigclip)
        bck_npix.append(len(np.where(clipdata.mask == False)[0]))
    bck_mean = np.array(bck_mean)
    bck_stdev = np.array(bck_stdev)
    bck_npix = np.array(bck_npix)
#
    phot = apphot(image, apertures)
    phot['bck_mean'] = bck_mean
    phot['bck_stdev'] = bck_stdev
    # Calculate bck-subtracted flux for each star plus photometry errors
    i = 0
    for aperture in apertures:
        area = aperture.area
        phot['flux_'+str(i)] = gain * phot['aperture_sum_'+str(i)] - bck_mean * area
        var = phot['flux_'+str(i)]/gain + area*bck_stdev**2 + area**2 * (bck_stdev)**2 / bck_npix
        phot['fluxerr_'+str(i)] = np.sqrt(var)
        i = i+1
    for col in phot.colnames:
        phot[col].info.format = '%.8g'
    #print(phot)
    phot.write(outfile, format='ascii', overwrite=True)


#--------------------------------------------------------------------
# Main script starts below
#--------------------------------------------------------------------

if len(sys.argv) < 2: 
    print(' Input FITS file name is required.')
    sys.exit()
    
infile = str(sys.argv[1])
if len(sys.argv) > 2:
    nsig = float(sys.argv[2])
else:
    nsig = 100.0
if len(sys.argv) > 3: 
    outtab = str(sys.argv[3])
else:
    outtab = 'simplephot.out'

# Divide science extension (in units of MJy/sr) by the number of MJy/sr producing 1 cps.
# Also get DQ array to enable exclusion of (badly) saturated sources
image = fits.getdata(infile, ext=1)
dqarr = fits.getdata(infile, ext=3)
print('Input image: ', infile)
fitsim = fits.open(infile)
scihead = fitsim[1].header
if 'PHOTMJSR' in scihead:
    photmjsr = fitsim[1].header['PHOTMJSR']
else:
    photmjsr = 1.0
image = image/photmjsr
# Define measurement radii here
myradii = [5.]
# Define pixel scale in arcsec here
mypixscale = 0.0656
# Define gain factor here
mygain = 1.6

source_list = find_sources(image, nsigma=nsig, roundness=1.0, sharpness=0.1)
phot_circ(image, dqarr, source_list, outtab, myradii, skip=5., w_ann=5., pixscale=mypixscale, gain=mygain)
