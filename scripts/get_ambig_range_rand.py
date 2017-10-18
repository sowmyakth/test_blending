"""Script computes points of separation of  ambiguously blended galaxies
from a random sample of 2 galaxy blends"""

import galsim
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from astropy.table import Table, Column
import lsst.afw.table
import lsst.afw.image
import lsst.afw.math
import lsst.meas.algorithms
import lsst.meas.base
import lsst.meas.deblender
import lsst.afw.detection


def make_table(num):
    names = ('NUM', 'flux1', 'flux2', 'hlr1', 'hlr2',
             'e1', 'e2', 'vmin', 'vmax', 'meas_flux')
    dtype = ['int'] + ['float'] * len(names) - 1
    cols = [range(num)] + [np.zeros(num)] * 4 + [np.zeros([num, 2])] * 5
    index_table = Table(cols, names=names, dtype=dtype)
    return index_table

def get_rand_param(table):

def get_range(hlr1, hlr2, low=0.9, high=1.8, num=30):
    """Returns distance between galaxies"""
    ratio = np.linspace(low, high, num)
    sum_val = (hlr1 + hlr2) / 0.2
    return sum_val * ratio


def get_hfpt(catalog, masked_image):
    """Returns heavy footprint of all objects in catalog"""
    hfpt = []
    for record in catalog:
        fpt = record.getFootprint()
        fpt.normalize()
        hfpt.append(lsst.afw.detection.makeHeavyFootprint(fpt, masked_image))
    return hfpt


def deblend_run(image_array, variance_array, psf_array):
    """Runs stack default settings on array and returns catalog and heavy
    footprint"""
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    config1 = lsst.meas.algorithms.SourceDetectionConfig()
    config1.tempLocalBackground.binSize = 128
    detect = lsst.meas.algorithms.SourceDetectionTask(schema=schema,
                                                      config=config1)
    deblend = lsst.meas.deblender.SourceDeblendTask(schema=schema)
    config1 = lsst.meas.base.SingleFrameMeasurementConfig()
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
    catalog2 = catalog.copy(deep=True)
    hfpt = get_hfpt(catalog2, masked_image)
    return catalog2, hfpt


def gal_im(hlr, flux, e, x0, y0, psf, nx, ny):
    """Return galsim image of galaxy with input parameters"""
    g = galsim.Convolve(galsim.Gaussian(flux=flux,
                                        half_light_radius=hlr).shear(e1=e[0],
                                                                     e2=e[1]),
                        psf)
    im = g.drawImage(scale=0.2, nx=nx, ny=ny,
                     offset=(x0, y0))
    return im


def get_var_arr(im, noise, snr):
    """Adds noise to image and returns array with variance value"""
    var = im.addNoiseSNR(noise, snr=500, preserve_flux=True)
    var_arr = var * np.ones_like(im.array)
    return var_arr


def get_sig_m(xx, yy, xy):
    """Return size sigma from input second moments"""
    det_Q = xx * yy - xy**2
    return np.abs(det_Q)**0.25


def met_rho(hfpt, hfpt_rest):
    """Returns blending parameter rho from input heavy footprint"""
    hfpt_rest.append(hfpt)
    denom = np.array([hfpt.dot(r) for r in hfpt_rest]).sum()
    return hfpt.dot(hfpt) / denom


def compute_val(x0s, hlr1, hlr2, flux1, flux2, e1, e2):
    """Computes blending parametrs for pair of blended input galaxy
    parametrs"""
    var = 1e-10
    nx, ny = 161, 161
    y0 = 0
    noise = galsim.GaussianNoise(sigma=var**0.5)
    psf = galsim.Gaussian(fwhm=0.7)
    psf_im = psf.drawImage(scale=0.2, nx=19, ny=19)
    psf_array = psf_im.array
    r = np.ones([len(x0s), 2]) * -1
    unit_dist = np.ones(len(x0s)) * -1
    num_detections = np.zeros(len(x0s))
    ret = []
    sig = [0, 0]
    for i, x0 in enumerate(x0s):
        im1_b_off = gal_im(hlr1, flux1, e1, 0, y0, psf, nx, ny)
        im2_b_off = gal_im(hlr2, flux2, e2, x0, y0, psf, nx, ny)
        im_b_on = im1_b_off + im2_b_off
        variance_array1 = get_var_arr(im1_b_off, noise, snr=500)
        variance_array2 = get_var_arr(im2_b_off, noise, snr=500)
        variance_array = get_var_arr(im_b_on, noise, snr=500)
        catalog1, hfpt1 = deblend_run(im1_b_off.array, variance_array1,
                                      psf_array)
        catalog2, hfpt2 = deblend_run(im2_b_off.array, variance_array2,
                                      psf_array)
        catalog, hfpt_blended = deblend_run(im_b_on.array, variance_array,
                                            psf_array)
        children = catalog[catalog['deblend_nChild'] == 0]
        if (len(catalog1) != 1) or (len(catalog2) != 1):
            print "multi detection in blending off: SKIPPED "
            continue
        elif len(children) == 0:
            print "No children in blending on"
            continue
        elif len(children) > 2:
            print "More than 2 children in blending on"
            continue
        elif np.isnan(children['base_SdssShape_flux']).any():
            print "Nan value found"
            continue
        tree = spatial.KDTree(zip([nx / 2, nx / 2 + x0],
                                  [ny / 2, ny / 2 + y0]))
        z = zip(children['base_SdssCentroid_x'],
                children['base_SdssCentroid_y'])
        match = tree.query(z, distance_upper_bound=10)
        select = ~np.isnan(match[0])
        if len(select) == 0:
            print "detected center do not match true object"
            continue
        # All good! Now procced with measurement
        num_detections[i] = len(children)
        ret.append(i)

        sig[0] = get_sig_m(catalog1['base_SdssShape_xx'],
                           catalog1['base_SdssShape_yy'],
                           catalog1['base_SdssShape_xy'])
        sig[1] = get_sig_m(catalog2['base_SdssShape_xx'],
                           catalog2['base_SdssShape_yy'],
                           catalog2['base_SdssShape_xy'])
        r[i][0] = met_rho(hfpt1[0], hfpt2)
        r[i][1] = met_rho(hfpt2[0], hfpt1)
        if len(children) == 1:
            sig0 = get_sig_m(children['base_SdssShape_xx'],
                             children['base_SdssShape_yy'],
                             children['base_SdssShape_xy'])
            unit_dist[i] = sig0 + sig[i][(match[1][select] + 1)%2]
        else:
            sigs = get_sig_m(children['base_SdssShape_xx'],
                             children['base_SdssShape_yy'],
                             children['base_SdssShape_xy'])
            unit_dist[i] = sigs.sum()
    s = np.array([x0s / unit_dist] * 2).T
    met = [r[ret], s[ret]]
    return met, num_detections[ret]


def get_transition_points(met, num_detections):
    """Returns transition point from 1 to 2 detections"""
    q1, = np.where(num_detections == 1)
    q2, = np.where(num_detections == 2)
    if (len(q1) == 0) or (len(q2) == 0):
        return [-1, -1], [-1, -1]
    min_vals = [min(met[i][q2].T[1]) for i in range(len(met))]
    max_vals = [max(met[i][q1].T[1]) for i in range(len(met))]
    return min_vals, max_vals


def main():
    num = 100
    tab = make_table(num)
    for i in range(num):
        x0s = get_range(tab['hlr1'][i], tab['hlr2'][i])
        met, det, flux = compute_val(x0s, tab['hlr1'][i], tab['hlr2'][i],
                                     tab['flux1'][i], tab['flux2'][i],
                                     tab['e1'][i], tab['e2'][i])
        v1, v2 = get_transition_points(met, det)
        tab['v1'][i] = v1
        tab['v2'][i] = v2
        tab['meas_flux'][i] = flux
