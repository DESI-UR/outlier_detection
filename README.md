# outlier_detection
Detection of outlier galaxy spectrum


Generating simulated galaxies:
### Generate DESI exposures in b, r, z cameras of various object types
1. generate_spectra.py
### Coadd observations from b, r, z into final spectrum and output the fits file
2. coadd_spectra.py (calls write_spectra.py to output final fits file)
### Inject SN templates into galaxies
3. inject_sne.py (TODO)
