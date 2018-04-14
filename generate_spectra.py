#################################################################
#                       generate_spectra.py                     #
# ------------------------------------------------------------- #
# Generate simulated DESI spectra for a variety of galaxy types #
# ------------------------------------------------------------- #
# Requires the DESI software:  https://github.com/desihub/      #
# Make sure you have the desi environment setup before running  # 
# ------------------------------------------------------------- #
# Author: Ryan Rubenzahl                                        #
# Last edit: 3/19/18                                            #
#################################################################

import os
import argparse

from astropy.io import fits
from astropy.table import Table

import desispec.io
import desisim.io
from desisim.obs import new_exposure
from desisim.scripts import quickgen
from desispec.scripts import group_spectra

p = argparse.ArgumentParser(description="Galaxy spectra generator with SN injection")
p.add_argument("-n", "--nspec", dest="nspec", default=100, type=int,
                help="Number of galaxies to simulate")
p.add_argument("-s", "--seed", dest="seed", default=555, type=int,
                help="Simulation seed")
p.add_argument("-f", "--flavor", dest="flavor", default='dark', action='store', type=str,
                choices=['dark', 'gray', 'grey', 'bright', 'bgs', 'mws', 'lrg', 
                            'elg', 'qso', 'std', 'arc', 'flat'], # all options
                help="Sky-brightness model and distribution of targets")
p.add_argument("-d", "--night", dest="night", default="20190101", type=str,
                help="Night of observation")
p.add_argument("-e", "--expid", dest="expid", default=0, type=int,
                help="Exposure ID number, can use to simulate more than one DESI exposure")
                    
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

print('')
print('Simulated raw data will be written to {}'.format(desisim.io.simdir()))
print('Pipeline will read raw data from {}'.format(desispec.io.rawdata_root()))
print('     (without knowing that it was simulated)')
print ('Pipeline will write processed data to {}'.format(desispec.io.specprod_root()))

# Generate the fibermap and specsim files needed by quickgen
# fibermap: table of simulated information about the position of each target in the DESI focal plane
# simspec:  "truth" table of spectra and the intrinsic properties of each object (e.g. redshift,
#               noiseless photometry, [OII] flux, etc.)
# If you want to use your own mix of objects, you need to write your own fibermap and specsim files
# following the particular format (rather than calling new_exposure)

# new_exposure generates random exposures of various types
output = new_exposure(args.flavor, nspec=args.nspec, seed=args.seed, night=args.night,
                        expid=args.expid, tileid=None, exptime=None)

# Check in on the data we just wrote
rawdata_dir = desispec.io.rawdata_root()
os.system('find {} | sort'.format(rawdata_dir))

fiberfile   = desispec.io.findfile('fibermap', night=args.night, expid=args.expid)
simspecfile = desisim.io.findfile('simspec',  night=args.night, expid=args.expid)

# Now we simulate the spectra using quickgen
quickgen.main(quickgen.parse(['--simspec', simspecfile,
                              '--fibermap', fiberfile])
              )


