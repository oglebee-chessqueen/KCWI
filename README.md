# KCWI: data manipulation codes to share between PC and linux
#
# Contains codes used to manipulate KCWI data cubes.
#
# Currently designed to look specifically at M57 data cubes:
#   - read in and store data, coordinates
#   - create narrow-band images around PNe-specific emission lines
#   - stitch together images taken of different regions of M57
#
# Codes created:
#   kcwi_datacube_general.py
#       - Was supposed to be a way to open and manipulate any KCWI data cubes, as long as data was reduced to the *_icube.fits level.
#       - Has been abandoned for now...
#   kcwi_image_shift.py
#       - Uses coordinates listed in FITS header (ra, dec) to match images to the same image.
#       - Had to also determine shifts of incorrect ra,dec listed in the images (due to mis-matching between guider and instrument during
#         commissioning)
#       - For M57, several different pointings were used over several nights, so this was a good way to tile the images
#         and match regions together -- very important for analysis later!
#   kcwi_scube_fix.py
#       - For nod+shuffle data cubes: During commissioning, n+s was tested but not optimized for sky field used.
#       - So, M57 n+s are contaminated with field stars. This program takes the *_scube.fits
#         files, finds where continuum sources are located, and replaces the elevated sky levels with the median per wavelength image.
#       - Then, I manually do *_ocube.fits - *_scube(fixed).fits and replace *_icube.fits with the result.
#   m57_kcwi_datacube_june.py
#       - Opens and reads in data and coordinates of all M57 data cubes.
#       - Includes lists of important emission lines.
#       - Creates select narrow-band images, ratio maps, and intensity per velocity bin.
#       - Good source of plotting code.
# 
# Codes in progress:
#   - Velocity/dynamics plots and analysis
#   - Electron density and temperature maps using line/ion ratios that the data covers
#
