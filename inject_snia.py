#################################################################
#                         inject_snia.py                        #
# ------------------------------------------------------------- #
# Inject SNe Ia templates into galaxy spectr                    #
# SNe Ia templates are from https://sne.space/download/         #
# ------------------------------------------------------------- #
# Requires the DESI software:  https://github.com/desihub/      #
# Make sure you have the desi environment setup before running  # 
# ------------------------------------------------------------- #
# Author: Ryan Rubenzahl                                        #
# Last edit: 4/04/18                                            #
#################################################################

import os
import json
import argparse

import numpy as np

from astropy.io import fits
from astropy.time import Time

from desispec.interpolation import resample_flux

p = argparse.ArgumentParser(description="Inject a supernova template into a galaxy spectrum")
p.add_argument("-g", "--galaxy", dest="args.galaxy_file", type=str,
                help="FITS file containing galaxy spectrum")
p.add_argument("-s", "--sntemp", dest="args.sn_template", type=str,
                help="FITS file containing supernova template")
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

# Read in supernova template
with open(args.sn_template, 'r') as f:
    sn_name = args.sn_template.split('.')[0].split('/')[-1]
    sn_data = json.load(f)
    
sn_type = sn_data[sn_name]['claimedtype'][0]['value']

# Grab a random epoch with sufficient observation
sn_waves  = []
sn_fluxes = []
sn_epochs = []

for sn_spectra in sn_data[sn_name]['spectra']:
    sn_flux_unit = sn_spectra['u_fluxes']
    sn_wave_unit = sn_spectra['u_wavelengths']
    sn_time_unit = sn_spectra['u_time']
    
    if sn_flux_unit == 'Uncalibrated':
        continue
        
    sn_time = sn_spectra['time']
    sn_spectra_data = np.array(sn_spectra['data'])
    sn_wave = np.array(sn_spectra_data[:,0], dtype=np.float32)
    sn_flux = np.array(sn_spectra_data[:,1], dtype=np.float32)
    
    if np.max(sn_wave) < 9000 or np.min(sn_wave) > 3000:
        continue

    sn_epoch = sn_peak_mjd - float(sn_time)
    sn_epochs.append(np.array(sn_epoch, dtype=np.float32))
    sn_waves.append(np.array(sn_wave, dtype=np.float32))
    sn_fluxes.append(np.array(sn_flux, dtype=np.float32))
    
rand_epoch = np.random.randint(len(sn_epochs))
sn_wave  = sn_waves[rand_epoch]
sn_flux  = sn_fluxes[rand_epoch]
sn_epoch = sn_epochs[rand_epoch]

# Get the peak overall flux of the supernova
sn_fluxes = np.array(sn_fluxes)
sn_epochs = np.array(sn_epochs)
sn_peak_spectra = np.concatenate(sn_fluxes[sn_epochs == np.min(np.abs(sn_epochs))])
sn_peak_flux = np.mean(sn_peak_spectra)

# Get lightcurve peak date
sn_peak_iso = sn_data[sn_name]['maxvisualdate'][0]['value']
sn_peak_mjd = Time(sn_peak_iso.replace('/', '-')).mjd

# Read in simulated galaxy spectra
gal_hdul = fits.open(args.galaxy_file)

gal_Z = gal_hdul[0].header['REDSHIFT']
gal_wave = gal_hdul[1].data
gal_flux = gal_hdul[2].data

gal_hdul.close()

# Redshift supernova template to the galaxy redshift
sn_wave *= (gal_Z + 1)

# Resample SN template to match galaxy spectra wavelength grid
sn_flux_ivar = None # TODO? add uncertainties
if sn_flux_ivar==None:
    sn_flux_resample = resample_flux(gal_wave, sn_wave, sn_flux, extrapolate=True)
    sn_flux_resample_ivar = None
else:
    sn_flux_resample, sn_flux_resample_ivar  = resample_flux(gal_wave, sn_wave, sn_flux, ivar=sn_flux_ivar)

# Scale supernova to galaxy
gal_flux_max = np.max(gal_flux)
sn_fluxratiorange = [0.1*gal_flux_max, gal_flux_max]

sn_fluxratio = np.random.uniform(sn_fluxratiorange[0], sn_fluxratiorange[1])
norm = sn_fluxratio / sn_peak_flux 

# Add together
sn_gal_flux = sn_flux_resample*norm + gal_flux

# Output to a new fits file
outfile = args.galaxy_file.split('.')[0] + "-with-SN" + sn_type + ".fits"

hdulist = fits.open(galaxy_file)

# Replace flux with supernova + galaxy flux
hdulist[2].data = sn_gal_flux

# Set header flag to indicate this spectra contains a supernova
hdulist[0].header['HAS_SN'] = sn_type

# Write to new fits file
hdulist.writeto(outfile)

hdulist.close()

print('written to {}'.format(outfile))
