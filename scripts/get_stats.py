"""Scrpit comapres true catalog to catalog of objects detected by DM and
identifies galaxies that are blended

To do
For now box for detected and truth is the same, but
edge effects may be presnt!!!
Shredded galaxies are not accounted for
Add fubction to get x y cord of true cat

The big galaxy!!!
check with noise

make getting image size a function
"""

import os
import numpy as np
from astropy.table import Table, Column
from scipy import spatial
from scipy.spatial.distance import cdist


def get_tue_cat(cat, Args):
    x1 = cat['dx'] / 0.2 + (4096 - 1) / 2.
    y1 = cat['dy'] / 0.2 + (4096 - 1) / 2.
    cond1 = (x1 >= int(Args.xmin) - 5) & (x1 < int(Args.xmax) + 5)
    cond2 = (y1 >= int(Args.ymin) - 5) & (y1 < int(Args.ymax) + 5)
    tru_cat = cat[np.where(cond1 & cond2)]
    return tru_cat


def get_det_cat(cat, Args):
    cond1 = (cat['deblend_nChild'] == 0)
    cond2 = (cat['base_GaussianCentroid_x'] >= int(Args.xmin))
    cond3 = (cat['base_GaussianCentroid_y'] >= int(Args.ymin))
    cond4 = (cat['base_GaussianCentroid_x'] < int(Args.xmax))
    cond5 = (cat['base_GaussianCentroid_y'] < int(Args.ymax))
    q, = np.where(cond1 & cond2 & cond3 & cond4 & cond5)
    det_cat = cat[q]
    return det_cat


def get_primary_detection(tru_cat, det_cat,
                          tolerance=5):
    x1 = tru_cat['dx'] / 0.2 + (4096 - 1) / 2.
    y1 = tru_cat['dy'] / 0.2 + (4096 - 1) / 2.
    z1_tree = spatial.KDTree(zip(x1, y1))
    x2 = det_cat['base_GaussianCentroid_x']
    y2 = det_cat['base_GaussianCentroid_y']
    z2 = zip(x2, y2)
    match = z1_tree.query(z2, distance_upper_bound=tolerance)
    col = Column(np.ones(len(x2)) * -1, 'Primary_Detection',
                 dtype='int')
    det_cat.add_column(col)
    col = Column(np.ones(len(x2)) * -1, 'Primary_Detection_gid',
                 dtype='int')
    det_cat.add_column(col)
    col = Column(np.ones(len(x2)) * -1, 'Prim_Det_Dist')
    det_cat.add_column(col)
    col = Column(np.zeros(len(tru_cat)), 'Detected',
                 dtype='int')
    tru_cat.add_column(col)
    tru_prim = match[1][~np.isinf(match[0])]
    dist_prim = match[0][~np.isinf(match[0])]
    det_cat['Primary_Detection'][~np.isinf(match[0])] = tru_prim
    det_cat['Primary_Detection_gid'][~np.isinf(match[0])] = tru_cat['grp_id'][tru_prim]
    det_cat['Prim_Det_Dist'][~np.isinf(match[0])] = dist_prim
    tru_cat['Detected'][tru_prim] = 1


def get_sig_m(xx, yy, xy):
    det_Q = xx * yy - xy**2
    return np.abs(det_Q)**0.25


def get_dist_ratio(tru_cat, det_cat):
    psf_sig = np.mean(get_sig_m(det_cat['base_SdssShape_psf_xx'],
                                det_cat['base_SdssShape_psf_yy'],
                                det_cat['base_SdssShape_psf_xy']))
    x1 = tru_cat['dx'][tru_cat['Detected'] != 1] / 0.2 + (4096 - 1) / 2.
    y1 = tru_cat['dy'][tru_cat['Detected'] != 1] / 0.2 + (4096 - 1) / 2.
    x2 = det_cat['base_GaussianCentroid_x']
    y2 = det_cat['base_GaussianCentroid_y']
    int_dist = cdist(zip(x2, y2), zip(x1, y1),
                     'euclidean')
    det_sigm = get_sig_m(det_cat['base_SdssShape_xx'],
                         det_cat['base_SdssShape_yy'],
                         det_cat['base_SdssShape_xy'])
    gal_sigm = np.sqrt(psf_sig**2 + (tru_cat['sigma_m'] / 0.2)**2)
    #t1 = np.resize(det_sigm**2, int_dist.shape)
    #t2 = np.resize(gal_sigm**2, int_dist.T.shape).T
    #new_unit_dist = np.sqrt(t1 + t2)
    t1 = np.resize(det_sigm, int_dist.shape)
    t2 = np.resize(gal_sigm, int_dist.T.shape).T
    new_unit_dist = t1 + t2
    return int_dist / new_unit_dist


def get_ambig(tru_cat, det_cat):
    ratio = get_dist_ratio(tru_cat, det_cat)
    ambig = np.where(ratio < 1)
    ambig_det = np.unique(ambig[0])
    ambig_tru = np.unique(ambig[1])
    col = Column(np.zeros(len(tru_cat)), 'ambig_blend', dtype='int')
    tru_cat.add_column(col)
    col = Column(np.zeros(len(det_cat)), 'ambig_blend', dtype='int')
    det_cat.add_column(col)
    no_det, = np.where(tru_cat['Detected'] != 1)
    tru_cat['ambig_blend'][no_det[ambig_tru]] = 1
    det_cat['ambig_blend'][ambig_det] = 1
    tru_cat['ambig_blend'][det_cat['Primary_Detection'][ambig_det]] = 1


def main(Args):
    # truth from catsim catalog (from wldebend)
    band = 'r'
    parentdir = os.path.abspath("..")
    tru_file = os.path.join(parentdir, 'data',
                            'wldeb_data/LSST_%s_trimmed.fits'%band)
    in_cat = Table.read(tru_file, format='fits',
                        hdu=1)
    tru_cat = get_tue_cat(in_cat, Args)
    # Results from DM stack
    det_file = os.path.join(parentdir, 'data',
                            'dm_output_%s.fits'%band)
    stack_cat = Table.read(det_file, format='fits',
                           hdu=1)
    det_cat = get_det_cat(stack_cat, Args)
    get_primary_detection(tru_cat, det_cat)
    get_ambig(tru_cat, det_cat)
    fname = os.path.join(parentdir, 'outputs',
                         'trucat_analyis_%s_%s.fit'%(band, Args.num))
    tru_cat.write(fname, format='fits', overwrite=True)
    fname = os.path.join(parentdir, 'outputs',
                         'detcat_analyis_%s_%s.fit'%(band, Args.num))
    det_cat.write(fname, format='fits', overwrite=True)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--xmin', default=0,
                        help="x cord of lower left pixel to run analysis[Default:0]")
    parser.add_argument('--ymin', default=0,
                        help="y cord of lower left pixel to run analysis[Default:0]")
    parser.add_argument('--xmax', default=4096,
                        help="x cord of top right pixel to run analysis[Default:100]")
    parser.add_argument('--ymax', default=4096,
                        help="y cord of top right pixel to run analysis[Default:100]")
    parser.add_argument('--num', default='0',
                        help="Number to denote output file name[Default:0]")
    args = parser.parse_args()
    main(args)
