"""Script to run LSST DM stack on a given fits image output from
WeakLensingDeblending package
To do:
install weaklensing deblending package in cori
make meas_args inputs read from blending package?
make psf more generic?
fix how the putput catalog nees to be saved
make variance array from noise read from image
check if output from debl package already has noise?

"""


import lsst.afw.table
import lsst.afw.image
import lsst.afw.math
import lsst.meas.algorithms
import lsst.meas.base
import lsst.meas.deblender
import lsst.meas.extensions.shapeHSM
import galsim
import math
import os
import numpy as np
from astropy.io import fits


def get_psf(Args):
    atmospheric_psf_fwhm = Args.zenith_psf_fwhm * Args.airmass**0.6
    if Args.atmospheric_psf_beta > 0:
        atmospheric_psf_model = galsim.Moffat(beta=Args.atmospheric_psf_beta,
                                              fwhm=atmospheric_psf_fwhm)
    else:
        atmospheric_psf_model = galsim.Kolmogorov(fwhm=atmospheric_psf_fwhm)
    lambda_over_diameter = 3600 * math.degrees(
                                               1e-10 * Args.central_wavelength / Args.mirror_diameter)
    area_ratio = Args.effective_area / (math.pi * (0.5 * Args.mirror_diameter)**2)
    obscuration_fraction = math.sqrt(1 - area_ratio)
    optical_psf_model = galsim.Airy(lam_over_diam=lambda_over_diameter,
                                    obscuration=obscuration_fraction)
    psf_model = galsim.Convolve(atmospheric_psf_model,
                                optical_psf_model)
    psf_size_pixels = 2 * int(math.ceil(10 * atmospheric_psf_fwhm / Args.pixel_scale))
    psf_image = galsim.Image(psf_size_pixels, psf_size_pixels,
                             scale=Args.pixel_scale)
    psf_model.drawImage(image=psf_image)
    return psf_image.array


def get_psf_file(psf_file):
    hdulist = fits.open(psf_file)
    psf_array = hdulist[0].data
    return psf_array


class meas_args(object):
    def __init__(self, band):
        if band == 'r':
            self.central_wavelength = 6199.52
            self.atmospheric_psf_beta = 0
            self.airmass = 1.2
            self.mirror_diameter = 8.36
            self.effective_area = 32.4
            self.pixel_scale = 0.2
            self.zenith_psf_fwhm = 0.7
        elif band == 'i':
            self.central_wavelength = 6199.52
            self.atmospheric_psf_beta = 0
            self.airmass = 1.2
            self.mirror_diameter = 8.36
            self.effective_area = 32.4
            self.pixel_scale = 0.2
            self.zenith_psf_fwhm = 0.67


def deblend_run(image_array, band,
                variance_array, psf_array, parentdir, Args):
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    config1 = lsst.meas.algorithms.SourceDetectionConfig()
    config1.tempLocalBackground.binSize = int(Args.binSize)
    # config1.doTempLocalBackground = True
    detect = lsst.meas.algorithms.SourceDetectionTask(schema=schema,
                                                      config=config1)
    deblend = lsst.meas.deblender.SourceDeblendTask(schema=schema)
    config1 = lsst.meas.base.SingleFrameMeasurementConfig()
    # config1.plugins.names.add('ext_shapeHSM_HsmShapeBj')
    # config1.plugins.names.add('ext_shapeHSM_HsmShapeLinear')
    # config1.plugins.names.add('ext_shapeHSM_HsmShapeKsb')
    config1.plugins.names.add('ext_shapeHSM_HsmShapeRegauss')
    config1.plugins.names.add('ext_shapeHSM_HsmSourceMoments')
    config1.plugins.names.add('ext_shapeHSM_HsmPsfMoments')
    measure = lsst.meas.base.SingleFrameMeasurementTask(schema=schema,
                                                        config=config1)
    image = lsst.afw.image.ImageF(image_array.astype(np.float32))
    variance = lsst.afw.image.ImageF(variance_array.astype(np.float32))
    masked_image = lsst.afw.image.MaskedImageF(image, None, variance)
    psf_image = lsst.afw.image.ImageD(psf_array.astype(np.float64))
    psf_kernel = lsst.afw.math.FixedKernel(psf_image)
    psf = lsst.meas.algorithms.KernelPsf(psf_kernel)
    exposure = lsst.afw.image.ExposureF(masked_image)
    exposure.setPsf(psf)
    # this is really just a factory for records, not a table
    table = lsst.afw.table.SourceTable.make(schema)
    detect_result = detect.run(table, exposure)
    # this is the actual catalog, but most of it's still empty
    catalog = detect_result.sources
    deblend.run(exposure, catalog)
    measure.run(catalog, exposure)
    catalog1 = catalog.copy(deep=True)
    det_file = os.path.join(parentdir, 'data',
                            'dm_output_%s_%sbin.fits'%(band, Args.binSize))
    catalog1.writeFits(det_file)


def main(Args):
    parentdir = os.path.abspath("..")
    band = 'i'  # 'r'
    im_file = os.path.join(parentdir, 'data',
                           'wldeb_data', 'LSST_%s_noise_trimmed.fits'%band)
    hdulist = fits.open(im_file)
    image_array = hdulist[0].data
    # args = meas_args(band)
    # psf_array = get_psf(args)
    # read PSF from file
    psf_file = os.path.join(parentdir, 'data',
                            'wldeb_data', 'mock_star.fits')
    psf_array = get_psf_file(psf_file)
    std = np.nanstd(image_array[0:50, 0:20])
    print "STD of noise", std
    std_data = 423.6466640862634  # STD from simualtion
    variance_array = np.ones_like(image_array) * std_data**2
    deblend_run(image_array, band,
                variance_array, psf_array, parentdir, Args)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--binSize', default=8,
                        help="tempLocalBackground binSize[Default:8]")
    args = parser.parse_args()
    main(args)
