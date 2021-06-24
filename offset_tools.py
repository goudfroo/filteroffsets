#!/usr/bin/env python
#
import sys
import os
import math
import numpy as np
from statistics import mode
from astropy import stats, table
from astropy.table import Table, Column, join_distance
from astropy.io import fits
from stsci.stimage import xyxymatch

def iterstat(image, nsigrej=5., maxiter=10, doprint=False,
              calcmode=False, lowlim='', highlim=''):
    """
    Perform iterative sigma-clipping statistics on image or array.
    Returns mean, sigma, median, minimum, maximum, and number of entries into a
    6-element array.
    Optionally also return mode as 7th parameter of output array.
    NOTE: recommend to use mode calculation only with Python >= 3.8.

    Parameters:
    -----------

    image:    Input image or array.
    nsigrej:  Optional parameter to designate "n" in "n times standard deviation". 
              Default is 5.
    maxiter:  Optional parameter to designate the max. number of clipping iterations. 
              Default is 10.
    doprint:  Optional parameter to print out results to STDOUT. Default is False.
    calcmode: Optional parameter to also calculate the mode. Default is False.
    lowlim:   Optional parameter to designate a lower limit to the values used in the 
              statistics. Default is no lower limit.
    highlim:  Optional parameter to designate a upper limit to the values used in the 
              statistics. Default is no upper limit.

    Paul Goudfrooij, June 2021

    """
    npimage = np.array(image)
    image_1d = npimage.flatten()
    if len(lowlim) > 0:
        llim = float(lowlim)
        ditchlow = np.where(image_1d > llim)
        image_1d = image_1d[ditchlow]
    if len(highlim) > 0:
        hlim = float(highlim)
        ditchhigh = np.where(image_1d < hlim)
        image_1d = image_1d[ditchhigh]
    clipdata = stats.sigma_clip(image_1d, sigma=nsigrej, maxiters=maxiter)
    mean, med, sig = stats.sigma_clipped_stats(image_1d, sigma=nsigrej, maxiters=maxiter)
    goods = np.where(clipdata.mask == False)[0]
    npix = len(goods)
    clipped = []
    for i in range(len(image_1d)):
        if i in goods:
            clipped.append(image_1d[i])
    minval = min(clipped)
    maxval = max(clipped)
    outarray = [mean, sig, med, minval, maxval, npix]
    if calcmode:
        modeval = mode(clipped)
        outarray.append(modeval)
        if doprint:
            print(' Statistics on '+str(image)+':')
            print(' Mean: {},  Sigma: {},  Median: {}'.format(mean,sig,med))
            print('    Min: {},  Max: {},  npix: {}'.format(minval,maxval,npix))
            print('    Mode: {}'.format(modeval))
            print(' ')
    else:
        if doprint:
            print(' Statistics on '+str(image)+':')
            print(' Mean: {},  Sigma: {},  Median: {}'.format(mean,sig,med))
            print('    Min: {},  Max: {},  npix: {}'.format(minval,maxval,npix))
            print(' ')
    return outarray


def titerstat(intab, col, rows=':', nsigrej=5., maxiter=10, doprint=False,
              calcmode=False, lowlim='', highlim=''):
    """
    Perform iterative sigma-clipping statistics on table column. 
    Returns mean, sigma, median, minimum, maximum, and number of entries into a
    6-element array.
    Optionally also return mode as 7th parameter of output array.
    NOTE: recommend to use mode calculation only with Python >= 3.8.

    Parameters:
    -----------

    intab:    Input table. ASCII, FITS, and several other formats are supported.
    col:      Column name of table. Use 'col' for ASCII tables without table headers.
    rows:     Optional parameter to subselect table rows. Use slice formatting 
              without brackets, like 'start:stop'. Default is all rows.
    nsigrej:  Optional parameter to designate "n" in "n times standard deviation". 
              Default is 5.
    maxiter:  Optional parameter to designate the max. number of clipping iterations. 
              Default is 10.
    doprint:  Optional parameter to print out results to STDOUT. Default is False.
    calcmode: Optional parameter to also calculate the mode. Default is False.
    lowlim:   Optional parameter to designate a lower limit to the values used in the 
              statistics. Default is no lower limit.
    highlim:  Optional parameter to designate a upper limit to the values used in the 
              statistics. Default is no upper limit.

    Paul Goudfrooij, June 2021

    """

    ttyp = os.popen('file '+intab)
    splittyp = ttyp.read().split()
    ftype = splittyp[1]
    if (ftype == 'ASCII'):
        t = Table.read(intab, format='ascii')
    elif (ftype == 'FITS'):
        t = Table.read(intab, format='fits')
    else:
        print('titerstat: automatically detecting table format. Hoping for the best.')
        t = Table.read(intab)
    incol = t[col]
    if len(rows) > 0:
        a = rows.split(":")[0]
        b = rows.split(":")[1]
        if len(a) > 0:
            index1 = int(a)
        else:
            index1 = 0
        if len(b) > 0:
            index2 = int(b)
        else:
            index2 = None
        incol = incol[index1:index2]
    if len(lowlim) > 0:
        llim = float(lowlim)
        ditchlow = np.where(incol > llim)
        incol = incol[ditchlow]
    if len(highlim) > 0:
        hlim = float(highlim)
        ditchhigh = np.where(incol < hlim)
        incol = incol[ditchhigh]
    clipdata = stats.sigma_clip(incol, sigma=nsigrej, maxiters=maxiter)
    mean, med, sig = stats.sigma_clipped_stats(incol, sigma=nsigrej, maxiters=maxiter)
    goods = np.where(clipdata.mask == False)[0]
    npix = len(goods)
    clipped = []
    for i in range(len(incol)):
        if i in goods:
            clipped.append(incol[i])
    minval = min(clipped)
    maxval = max(clipped)
    outarray = [mean, sig, med, minval, maxval, npix]
    if calcmode:
        modeval = mode(clipped)
        outarray.append(modeval)
        if doprint:
            print(' Statistics on column '+col+':')
            print(' Mean: {},  Sigma: {},  Median: {}'.format(mean,sig,med))
            print('    Min: {},  Max: {},  npix: {}'.format(minval,maxval,npix))
            print('    Mode: {}'.format(modeval))
            print(' ')
    else:
        if doprint:
            print(' Statistics on column '+col+':')
            print(' Mean: {},  Sigma: {},  Median: {}'.format(mean,sig,med))
            print('    Min: {},  Max: {},  npix: {}'.format(minval,maxval,npix))
            print(' ')
    return outarray


def coordmatch(tab1, tab2, outtab, matchcols1, matchcols2, maxdist,
               minsep=0.0, incols1='', incols2=''):

    intab1 = str(tab1)
    intab2 = str(tab2)
    outfile = str(outtab)
    distmax = float(maxdist)

    if len(matchcols1) > 2:
        print(' Need 1 or 2 columns for the matching. Exiting coordmatch.')
        sys.exit()
    if len(matchcols1) != len(matchcols2):
        print(' Mismatch in the number of columns to match ({} vs. {}).'.format(len(matchcols1),len(matchcols2)))
        print(' Exiting coordmatch.')
        sys.exit()
          
    t1 = Table.read(intab1, format='ascii')
    t2 = Table.read(intab2, format='ascii')

# Default in this case is not to include any columns from the input tables
# (other than the coordinate columns that is). So commenting this out:
#if "incols1" not in locals():
#    incols1 = []
#    for col in t1.colnames:
#          if col not in matchcols1:
#             incols1.append(col)
#if "incols2" not in locals():
#    incols2 = []
#    for col in t2.colnames:
#          if col not in matchcols2:
#             incols2.append(col)

    if (len(matchcols1) == 1): 
              x1 = t1[matchcols1[0]]
              x2 = t2[matchcols2[0]]
              tt1 = Table([x1], names=['coord'])
              tt2 = Table([x2], names=['coord'])
              t12 = table.join(tt1, tt2, join_funcs={'coord': join_distance(distmax)})
              dist = abs(t12['coord_1']-t12['coord_2'])
              # Then subselect only matches more distant than "separation"
              distsel = np.where(dist > minsep)
              thedist = dist[distsel]
              thecoordid = coord_id[distsel]
              thecoord1 = coord_1[distsel]
              thecoord2 = coord_2[distsel]
              thetab = Table([thecoordid, thecoord1, thecoord2],
                              names=['coord_id', matchcols1[0]+'_1',
                                     matchcols2[0]+'_2'])
              # Then search for those occurrences to get the input table rows
              if ((len(incols1) > 0) or (len(incols2) > 0)):
                    for i in range(len(thecoord1)):
                        rows1 = []
                        rows2 = []
                        for j in range(len(x1)):
                              if (x1[j] == thecoord1[i]):
                                  rows1.append(j)
                        for k in range(len(x2)):
                              if (x2[k] == thecoord2[i]):
                                  rows2.append(k)
                    # Then get the desired rows/columns from the input tables
                    if (len(incols1) > 0):
                        for i in range(len(incols1)):
                              oldcol = incols1[i]
                              newcol = []
                              for j in rows1:
                                  newcol.append(t1[oldcol][j])
                              thetab.add_column(newcol, name=oldcol+'_1')
                    if (len(incols2) > 0):
                        for i in range(len(incols2)):
                              oldcol = incols2[i]
                              newcol = []
                              for j in rows2:
                                  newcol.append(t2[oldcol][j])
                              thetab.add_column(newcol, name=oldcol+'_2')
              thetab['coord_id'].info.format = '.4d'
              thetab.write(outtab, format='ascii', overwrite=True)  
    else:
              x1 = t1[matchcols1[0]]
              y1 = t1[matchcols1[1]]
              x2 = t2[matchcols2[0]]
              y2 = t2[matchcols2[1]]
              xy1 = np.stack([x1, y1]).T
              xy2 = np.stack([x2, y2]).T
              matches = xyxymatch(xy2, xy1, tolerance=distmax, separation=minsep)
              matchtab = Table(matches)
              dx = matchtab['input_x']-matchtab['ref_x']
              dy = matchtab['input_y']-matchtab['ref_y']
              dist = np.sqrt(dx**2 + dy**2)
              matchtab.add_column(dx, name='dx')
              matchtab.add_column(dy, name='dy')
              matchtab.add_column(dist, name='distance')
              # Now add the selected rows from the desired columns of input table #1
              # (in this case this is likely irrelevant, but still keeping it in here)
              if (len(incols1) > 0):
                    goodid1 = matchtab['input_idx']
                    for i in range(len(incols1)):
                        oldcol = incols1[i]
                        newcol = []
                        for j in goodid1:
                              newcol.append(t2[oldcol][j])
                        matchtab.add_column(newcol, name='in_'+oldcol)
                        matchtab['in_'+oldcol].info.format = '.4g'
              # Now add the selected rows from the desired columns of ref. table #2
              # (in this case this is likely irrelevant, but still keeping it in here)
              if (len(incols2) > 0):
                    goodid2 = matchtab['ref_idx']
                    for i in range(len(incols2)):
                        oldcol = incols2[i]
                        newcol = []
                        for j in goodid2:
                              newcol.append(t1[oldcol][j])
                        matchtab.add_column(newcol, name='ref_'+oldcol)
                        matchtab['ref_'+oldcol].info.format = '.4g'
              # Enter column formats and write output table
              matchtab['dx'].info.format = '8.4f'
              matchtab['dy'].info.format = '8.4f'
              matchtab['distance'].info.format = '9.4f'
              matchtab.write(outtab, format='ascii', overwrite=True)
