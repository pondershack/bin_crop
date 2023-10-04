
"""
This program creates a copy of an existing calibration image file in order to work around an image format incompatibility,\
    making it usable by EXOTIC and AstroImageJ. 
    It crops and bins the images, and adjusts FITS keywords where needed, using the fits module from astropy.
    The new calibration image is in the 32-bit floating point format.
"""

import sys
import getopt
from astropy.io import fits
import astropy.io
import numpy as np
from typing import List


def get_images():

    """
    This function 
    * uses getopt to parse inputs into three variables
    * gives an error message when needed.
    """

    global calibration_image, science_image, new_calibration_image

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "i:t:o:", ["calibration_image=", "science_image=", "new_calibration_image="])
    except getopt.GetoptError as err:
        print(err)
        opts = []
        #exit()
    
    for opt, arg in opts:
        if opt in ['-i', '--calibration_image']:
            calibration_image = arg
        elif opt in ['-t', '--science_image']:
            science_image = arg
        elif opt in ['-o', '--new_calibration_image']:
            new_calibration_image = arg

def rebin(arr, new_shape):

            """
            This function reshapes array into a higher dimension and then takes the means over the appropriate new axes
            Bins larger array into smaller one assuming that each dimension of the smaller is a factor
            of the corresponding dimension in the larger.
            """

            shape = (new_shape[1], arr.shape[0] // new_shape[1], new_shape[0], arr.shape[1] // new_shape[0])
            return arr.reshape(shape).mean(-1).mean(1)



class FITSImage:

    def __init__(self):
        self.name : str = None
        self.imageHDU : int = 0

    def find_image_hdu(self):
        with fits.open(str(self.name), mode='readonly') as list:
            """
            This routine finds the image hdu
            """
            i = -1

            while True:
                try:
                    i += 1
                    value = list[i].header['EXPOSURE']
                    self.imageHDU = i
                except KeyError:
                    continue
                except IndexError:
                    # went through all hdus
                    break

    def view_file(self):
        with fits.open(str(self.name), mode='readonly') as list:
            """
            This routine just views the image
            """
            list.info()
            i = -1

            while True:
                try:
                    i += 1
                    hdr = list[i].header
                    print(repr(hdr))
                except IndexError:
                    break

    def copy_file(self):
        with fits.open(str(self.name), mode='readonly') as list:
            """
            This routine
            * copies the image to new_calibration_image
            """

            list.writeto(str(new_calibration_image), overwrite=False)
            #creates a new calibration file which is an exact copy of the old one
            #overwriting files is turned off for safety, this might cause an error, but it's simple to fix

    def float32_convert(self):
        with fits.open(str(self.name), mode='update') as list:
            """
            This routine
            * converts image to 32 float permanently
            """

            data = list[self.imageHDU].data
            #opens array of image

            list[self.imageHDU].data = data.astype(np.float32)
            #converts to float 32 
            #Note: when editing data or hdr0, use calibration_list[1].data or calibration_list[0].header instead

            list.flush()
            #changes the original file

    def get_keywords(self):
        with fits.open(str(self.name), mode='readonly') as list:
            """
            This routine
            * stores dimensions of image
            * stores values of FRAMEX, FRAMEY, BINNING, DATAMAX for image
            * prints errors if needed
            """

            data = list[self.imageHDU].data

            Y, X = data.shape
            #gets dimensions

            hdr = list[self.imageHDU].header

            try:
                FRAMEX, FRAMEY, BINNING, DATAMAX = hdr['FRAMEX'], hdr['FRAMEY'], hdr['BINNING'], hdr['DATAMAX']
            except TypeError as err:
                #TypeError is the correct choice when values are missing
                print(err)
                FRAMEX, FRAMEY, BINNING, DATAMAX = None, None, None, None
        
        return [X, Y, FRAMEX, FRAMEY, BINNING, DATAMAX]
    
    def left_crop(self):
        with fits.open(str(new_calibration_image), mode='update') as new_calibration_list:
            """
            This routine:
            * crops the image horizontally
            * adjusts FRAMEX keyword
            * saves changes to new_calibration_image
            * only the new_calibration_image is cropped, so no need for self.name
            """

            data = new_calibration_list[self.imageHDU].data

            new_calibration_list[self.imageHDU].data = new_calibration_list[self.imageHDU].data[0:calibrationY, ((scienceFRAMEX - calibrationFRAMEX) // calibrationBINNING) : (calibrationX-((scienceFRAMEX - calibrationFRAMEX) // calibrationBINNING))]
            #crops image

            hdr = new_calibration_list[self.imageHDU].header
            hdr['FRAMEX'] = (scienceFRAMEX - calibrationFRAMEX)
            # changes the value of the FRAMEX keyword to the necessary value

            new_calibration_list.flush()

    def binning(self):
        with fits.open(str(new_calibration_image), mode='update', memmap = False) as new_calibration_list:
            # memmap = False prevents weird errors
            """
            This whole section:
            * reshapes new_calibration_image slightly to make binning possible
            * bins the new_calibration_image to that of the science_image
            * adjusts BINNING and DATAMAX keywords
            * saves changes to new_calibration_image
            * binning is done on new_calibration_image so self.name is unnecessary
            """

            new_calibration_list[self.imageHDU].data = new_calibration_list[self.imageHDU].data[0:(new_calibrationY - (new_calibrationY % scienceY)), 
                                                                        0:(new_calibrationX - (new_calibrationX % scienceX))]

            """
            This code makes it so that each dimension in the science_image array is a factor of the corresponding dimension
            in the new_calibration_image array
            """

            data = new_calibration_list[self.imageHDU].data

            new_calibration_list[self.imageHDU].data = rebin(data, (scienceX, scienceY))
            # does the binning

            hdr = new_calibration_list[self.imageHDU].header
            hdr['BINNING'] = (scienceBINNING)
            hdr['DATAMAX'] *= (scienceBINNING ** 2)
            # changes the value of the BINNING and DATAMAX keyowrds to the necessary value
   

if __name__ == '__main__':
    
    calibration_image = None
    science_image = None
    new_calibration_image = None
    # create variables for image names as strings

    get_images()

    print('calibration_image: {}'.format(calibration_image))
    print('science_image: {}'.format(science_image))
    print('new_calibration_image: {}'.format(new_calibration_image))

    image = FITSImage()
    image.name = calibration_image
    image.find_image_hdu()
    calibrationX, calibrationY, calibrationFRAMEX, calibrationFRAMEY, calibrationBINNING, calibrationDATAMAX = image.get_keywords()
    image.copy_file()

    image.name = science_image
    image.find_image_hdu()
    scienceX, scienceY, scienceFRAMEX, scienceFRAMEY, scienceBINNING, scienceDATAMAX = image.get_keywords()

    image.name = new_calibration_image
    image.find_image_hdu()
    image.float32_convert()

    if calibrationFRAMEX != scienceFRAMEX:
        leftCrop = True
        # these variables are used for the final message at the end
        image.left_crop()

    elif calibrationFRAMEX == scienceFRAMEX:
        leftCrop = False

    if calibrationFRAMEY != scienceFRAMEY:
        print("Top needs to be cropped")
        topCrop = True
        sys.exit(1)

    elif calibrationFRAMEY == scienceFRAMEY:
        topCrop = False

    if calibrationBINNING == scienceBINNING:
        print("No binning is needed")
        Binning = 0

    elif calibrationBINNING == 1:
        Binning = 1
        new_calibrationX, new_calibrationY, new_calibrationFRAMEX, new_calibrationFRAMEY, new_calibrationBINNING, new_calibrationDATAMAX = image.get_keywords()
        image.binning()

    else:
        print('Special binning is needed')
        Binning = 2

    new_calibrationX, new_calibrationY, new_calibrationFRAMEX, new_calibrationFRAMEY, new_calibrationBINNING, new_calibrationDATAMAX = image.get_keywords()
    image.view_file()

    print("This program created a copy of an existing calibration image file in order to work around an image format incompatibility,\
    making it usable by EXOTIC and AstroImageJ. The new calibration image is in the 32-bit floating point format.")

    if leftCrop == True:
        print("This program cropped the left and right edges of the calibration image to match those of the pre-binned science image,\
        adjusting the necessary keywords in the proper HDU. Final cropping value:{}".format(new_calibrationFRAMEX))
    elif leftCrop == False:
        print("This program did not crop the left or the right")

    if topCrop == True:
        print("This program cropped the top and bottom edges of the calibration image to match those of the pre-binned science image,\
        adjusting the necessary keywords in the proper HDU. Final cropping value:{}".format(new_calibrationFRAMEY))
    elif topCrop == False: 
        print("This program did not crop the top or the bottom")

    if Binning == 1:
        print("This program binned the calibration image to match the science image. Final binning value:{}".format(new_calibrationBINNING))

    elif Binning == 0:
        print("No binning was executed")

    if Binning == 2:
        print("This program binned the calibration image to match the science image. Final binning value:{}".format(new_calibrationBINNING))

    """
    Comments:
    * What about the keyword for degrees of sky per pixel?
    * Errors need testing
    * program should work for any left cropping, testing is needed
    * doesn't crop vertically
    * should be able to bin anything, as long as calibration binning is 1, testing is needed
    * definitely is able to bin into 3 by 3
    """

    #  calling the program example:
    #  python bin_crop.py -i _dark60_16.fits -t _image300.fits -o NEWFILE.fits