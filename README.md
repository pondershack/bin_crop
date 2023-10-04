This program creates a copy of an existing calibration image file in order to work around an image format incompatibility,\
    making it usable by EXOTIC and AstroImageJ. 
    It crops and bins the images, and adjusts FITS keywords where needed, using the fits module from astropy.
    The new calibration image is in the 32-bit floating point format.
