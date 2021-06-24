#!/usr/bin/env python
#
# Syntax: python filteroffsets.py <intab1> <intab2> <outtab> <maxdist> 
#                nsigrej=5. matchcol1='xcenter' matchcol2='ycenter' 
#                minsep=0.0 incols1='' incols2=''
#
# Matches table rows between two input tables to within a distance of <maxdist>.
# The X and Y coordinates of the matched sources are put into output table
# <outtab>. The script then calculates kappa-sigma clipping statistics on the
# offsets in X and Y, which are columns "dx" and "dy" in table <outtab>. 
# Results of the statistics are printed to STDOUT, and they also go into array
# variables "dxstats" and "dystats" when used within a python shell or ipython.
#
# Note that the second input table <intab2> is treated as "reference" table.
#
# Other (optional) parameters:
# ----------------------------
#
#   nsigrej: number of sigmas to be used in clipping (defaulted to 5.)
# matchcol1: Name of column with X coordinates in input tables
#            (defaulted to 'xcenter')
# matchcol2: Name of column with y coordinates in input tables
#            (defaulted to 'ycenter')
#    minsep: Minimum separation between matching sources
#            (defaulted to 0.0)
#   incols1: Columns from table <intab1> to be copied to <outtab>
#            (defaulted to no columns in addition to the coordinates)
#   incols2: Columns from table <intab2> to be copied to <outtab>
#            (defaulted to no columns in addition to the coordinates)
#
# P. Goudfrooij
# last update: 23 June 2021
#
import sys
import os
import math
import numpy as np
from astropy import units, stats, table
from astropy.io import fits
from astropy.table import Table, Column
from offset_tools import coordmatch, titerstat

#------------------- Input parameters -------------------------
if len(sys.argv) < 5: 
    print('Syntax:')
    print('python filteroffsets.py <intab1> <intab2> <outtab> <maxdist>')
    print('         nsigrej=5. matchcol1="xcenter" matchcol2="ycenter" minsep=0.0')
    print('         cols1="" cols2=""')
    print('Note: Input table #2 is used as reference table.')
    sys.exit()
    
intab1 = str(sys.argv[1])
intab2 = str(sys.argv[2])
outtab = str(sys.argv[3])
maxdist = float(sys.argv[4])
if len(sys.argv) > 5:
    nsig = float(sys.argv[5])
else:
    nsig = 5.0
if len(sys.argv) > 6:
    matchcol1 = str(sys.argv[6])
else:
    matchcol1 = 'xcenter'
if len(sys.argv) > 7:
    matchcol2 = str(sys.argv[7])
else:
    matchcol2 = 'ycenter'
if len(sys.argv) > 8:
    mindist = float(sys.argv[8])
else:
    mindist = 0.0
if len(sys.argv) > 9:
    cols1 = sys.argv[9]
else:
    cols1 = ''
if len(sys.argv) > 10: 
    cols2 = sys.argv[10]
else:
    cols2 = ''

# Concatanate the two coordinate columns into one variable
matchcols = [matchcol1,matchcol2]

# Coordinates matched between the two input tables go into output table
coordmatch(intab1, intab2, outtab, matchcols, matchcols, maxdist,
           minsep=mindist, incols1=cols1, incols2=cols2)
# Do iterative statistics on x and y offsets between two sets of coordinates.
# Results go into output variables "dxstats" and "dystats".
dxstats = titerstat(outtab, 'dx', nsigrej=nsig, doprint=True)
dystats = titerstat(outtab, 'dy', nsigrej=nsig, doprint=True)

