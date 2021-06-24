# JWST Filter Offsets Workflow
This repository contains scripts that can be used for determining **Filter-to-Filter Offsets** from a set of JWST images that were taken with different filters without any spatial dithers executed between the individual images.

## Installation
This package has some external dependencies. Specifically, it requires the packages NumPy, AstroPy, statistics, and stsci.stimage. Python version 3.8 or later is preferred but not required if one selects the default settings of the scripts in this package. This package will work fine within the conda environment "jwst_fpa" (used for commissioning program NIS-013, see
[this repo](https://github.com/tonysohn/jwst_fpa)) or the general
[niriss-commissioning repo](https://github.com/spacetelescope/niriss-commissioning). The instructions below assume that one of those environments have been installed on your machine and that you have loaded into that environment.

# Usage
There are two main scripts and one supplemental script (`offset_tools.py`, which includes utility functions). Below are brief descriptions on the use of the two main scripts.

(1) `simplephot.py` - This script takes a `_cal.fits` (or `_rate.fits`) file as input, then establishes a catalog of point-like sources in the image, then filters out the sources that either contain more than one saturated pixel or are located in a region with bad "do not use" pixels, and performs simple astrometry and aperture photometry on the remaining sources. The script is run as follows:

```
python simplephot.py <input FITS file> <nsigma> <output table>
```

here, the optional parameter <nsigma> denotes the N in "N * sigma" for the source detection threshold (sigma being the standard deviation of the sky background in the image), which is defaulted to 100. The optional parameter <output table>, which is defaulted to "simplephot.out", denotes the name of the (ASCII-format) output table that holds the astrometry and photometry values.

(2) `filteroffsets.py` - This script takes two tables coming from simplephot.py as input (typically one for each of two passbands. the second input table is used as "reference" table, and this table should represent the passband that is to be used as "zero offset" reference), then matches source coordinates from those two tables, and puts the X and Y coordinates of the matched sources in an output table. Finally, it performs iterative kappa-sigma clipping statistics on the differences of the X and Y coordinates from the two tables and outputs the results to the terminal. The script is run as follows:

```
python filteroffsets.py <intab1> <intab2> <outtab> <maxdist>
```

where <intab1> and <intab2> are the two input tables (which are envisaged to be output tables from `simplephot.py`), <outtab> is the output table that will contain the X and Y coordinates of the sources matched between the two input tables, and <maxdist> is the maximum distance between sources in the two input tables to be included as a "matched" source. Note that there are 6 additional optional input parameters to the `filteroffsets.py` script, which are defaulted to values deemed appropriate for this commissioning activity. Information on these optional parameters can be found in the script itself.


## Contributing

Any changes to code should be made on a feature branch. The workflow should be:

1. Create a branch
2. Commit changes on branch
3. Submit pull request for branch to be merged into master

Detailed instructions can be found in the 
[contribution guidelines](https://github.com/spacetelescope/niriss-commissioning/blob/master/CONTRIBUTING.md).

## Code of Conduct

All users are expected to adopt a 
[code of conduct](https://github.com/spacetelescope/niriss-commissioning/blob/master/CODE_OF_CONDUCT.md) 
that ensures a productive, respectful environment for all contributors and 
participants. Any issues or violations should be reported to conduct@stsci.edu. 

## Questions or Issues

Any code problems or questions should be noted by 
[opening an issue](https://github.com/spacetelescope/niriss-commissioning/issues).

If you have git questions, see 
[these resources](https://github.com/spacetelescope/niriss-commissioning/blob/master/CONTRIBUTING.md#Resources).


