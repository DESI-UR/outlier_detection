#################################################################
#                        coadd_spectra.py                       #
# ------------------------------------------------------------- #
# Coadd simulated DESI spectra: combine exposures from multiple #
# wavelengths (cameras) into one overall spectrum               #
# ------------------------------------------------------------- #
# Requires the DESI software:  https://github.com/desihub/      #
# Make sure you have the desi environment setup before running  # 
# ------------------------------------------------------------- #
# Author: Ryan Rubenzahl                                        #
# Last edit: 3/19/18                                            #
#################################################################

import os
import argparse

import numpy as np

from astropy.io import fits
from astropy.table import Table

import desispec.io
import desisim.io
from desispec.coaddition import Spectrum
from desispec.resolution import Resolution

from write_spectra import write_coadd_spectra

p = argparse.ArgumentParser(description="Galaxy spectra generator with SN injection")
p.add_argument("-o", "--outdir", dest="outdir", default="./output/", type=str,
                help="Path to output directory")
p.add_argument("-f", "--flavor", dest="flavor", default='dark', action='store', type=str,
                choices=['dark', 'gray', 'grey', 'bright', 'bgs', 'mws', 'lrg', 
                            'elg', 'qso', 'std', 'arc', 'flat'], # all options
                help="Sky-brightness model and distribution of targets")
p.add_argument("-d", "--night", dest="night", default="20190101", type=str,
                help="Night of observation")
p.add_argument("-e", "--expid", dest="expid", default=0, type=int,
                help="Exposure ID number, can use to simulate more than one DESI exposure")
p.add_argument("-n", "--nadd", dest="nadd", default=None, type=int,
                help="Number of spectra to add (if None will default to all)")
args = p.parse_args()

def check_env():
	for env in ('DESIMODEL', 'DESI_ROOT', 'DESI_SPECTRO_SIM', 'DESI_SPECTRO_DATA',
								'DESI_SPECTRO_REDUX', 'SPECPROD', 'PIXPROD'):
		if env in os.environ:
			print('{} environment set to {}'.format(env, os.getenv(env)))
		else:
			print('Required environment variable {} not set!'.format(env))

# Make sure environments are setup properly
check_env()

# Reassign output directories so we don't mess anything up in the standard DESI integration test
os.environ["SPECPROD"] = "example"
os.environ["PIXPROD"]  = "example"
rawdata_dir = desisim.io.simdir()
os.environ["DESI_SPECTRO_DATA"] = rawdata_dir

print('Simulated raw data will be read from {}'.format(desisim.io.simdir()))
print('     (without knowing that it was simulated)')

# Read in fibermap and simspec files
fiberfile   = desispec.io.findfile('fibermap', night=args.night, expid=args.expid)
simspecfile = desisim.io.findfile('simspec',  night=args.night, expid=args.expid)

print ('Reading simspec file {}.'.format(simspecfile))
hdu = fits.open(simspecfile)
meta = Table(hdu['TRUTH'].data)
hdu.close()

print('Reading fibermap file {}.'.format(fiberfile))
hdu = fits.open(fiberfile)
hdu.info()
fibermap = Table(hdu['FIBERMAP'].data)
hdu.close()

nspec = len(fibermap)

# Now read in simulated spectra
cameras = ['b', 'r', 'z']
cframes = {camera: [] for camera in cameras}

# read the data for each camera and store in dictionary
for cam in cameras:
    for i in range(0, nspec, 500):
        camera = cam +  "{}".format(i/500)

        cframefile = desispec.io.findfile('cframe', night=args.night, expid=args.expid, camera=camera)
        print('Reading {}'.format(cframefile))
        try:
            cframes[cam] = desispec.io.frame.read_frame(cframefile)
        except IOError:
            print('cframefile not found {}'.format(cframefile))
            continue

# How many spectra to add
if args.nadd == None:
    nadd = nspec
elif args.nadd > nspec:
    nadd = nspec
else:
    nadd = args.nadd

# Now, coadd the b, r, and z spectra into one final spectrum
for n in range(nadd):
    # First, get the overall wavelength range of the spectrum
    # Individual cameras have 0.5 A resolutions
    # We make the final coadded spectrum uniformly spaced with 1.0 A bins
    global_wavelength_grid = np.arange(cframes['b'].wave[0],
                                       cframes['z'].wave[-1] + 0.5,
                                       1.0, dtype=np.float32)
    
    coadd_all_bands = Spectrum(global_wavelength_grid)
    for band in cframes:
        band_specobj = Spectrum(cframes[band].wave,
                                flux=cframes[band].flux[n],
                                ivar=cframes[band].ivar[n],
                                resolution=Resolution(cframes[band].resolution_data[n]))

        coadd_all_bands += band_specobj 

    # Finalize the overall spectrum
    coadd_all_bands.finalize()
    print('{}/{} coadded'.format(n+1, nadd))

    # And write the output
    outfile = args.outdir + "spectra-{:05d}.fits".format(n)
    write_coadd_spectra(outfile, coadd_all_bands, fibermap[n], simspecfile)
    print('written to {}'.format(outfile))
