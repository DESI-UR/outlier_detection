#################################################################
#                        write_spectra.py                       #
# ------------------------------------------------------------- #
# Write the coadded spectra to an output fits file              #
# ------------------------------------------------------------- #
# Requires the DESI software:  https://github.com/desihub/      #
# Make sure you have the desi environment setup before running  # 
# ------------------------------------------------------------- #
# Author: Ryan Rubenzahl                                        #
# Last edit: 3/19/18                                            #
#################################################################
import os

from astropy.io import fits

from desispec.io.util import fitsheader

def write_coadd_spectra(outfile, coadd_spectra, fiberdata, metadata, simspecfile, units=None):
    """
    Write the output of the full brz coadded spectra to a fits file.
    Include the truth and metadata from fiberfile and simspecfile
    """
    
    # Make parent directory if necessary
    dir, base = os.path.split(outfile)
    if not os.path.exists(dir):
        os.makedirs(dir)
        
    # Create HDUs from the data
    all_hdus = fits.HDUList()
    
    sshdu = fits.open(simspecfile)
    
    metadict = {'NIGHT': sshdu[0].header['NIGHT'],
                'EXPID': sshdu[0].header['EXPID'],
                'TILEID': sshdu[0].header['TILEID'],
                'HAS_SN': False}
   
    # Add other metadata
    # TODO: dump all metadata from Table into metadict
    metadict['REDSHIFT'] = metadata['REDSHIFT']
    metadict['OBJTYPE']  = metadata['OBJTYPE']

    # metadata goes in empty primary HD
    hdr = fitsheader(metadict)
    all_hdus.append(fits.PrimaryHDU(header=hdr))

    # Wavelength data
    hdu = fits.ImageHDU(name='WAVE')
    hdu.header["BUNIT"] = "Angstrom"
    hdu.data = coadd_spectra.wave.astype("f8")
    all_hdus.append(hdu)
    
    # Flux data
    hdu = fits.ImageHDU(name='FLUX')
    if units is None:
        hdu.header["BUNIT"] = "1e-17 erg/(s cm2 Angstrom)"
    else:
        hdu.header["BUNIT"] = units
    hdu.data = coadd_spectra.flux.astype("f4")
    all_hdus.append(hdu)
    
    # Variance data
    hdu = fits.ImageHDU(name="IVAR")
    hdu.data = coadd_spectra.ivar.astype("f4")
    all_hdus.append(hdu)
    
    # Mask
    if coadd_spectra.mask is not None:
        hdu = fits.CompImageHDU(name="MASK")
        hdu.data = coadd_spectra.mask.astype(np.uint32)
        all_hdus.append(hdu)
        
    # Resolution data
    if False:#coadd_spectra.resolution is not None:
        hdu = fits.ImageHDU(name="RESOLUTION")
        hdu.data = coadd_spectra.resolution.astype("f4")
        all_hdus.append(hdu)
        
    try:
        all_hdus.writeto("{}.tmp".format(outfile), 
                         overwrite=True, checksum=True)
    except TypeError:
        all_hdus.writeto("{}.tmp".format(outfile), 
                         clobber=True, checksum=True)
    os.rename("{}.tmp".format(outfile), outfile)
    
    return outfile
