import os
import numpy as np
from astropy.table import Table, vstack
"""Script creates pair of galaxies in the format of catsim galaxy catalog.
The central galaxy is picked at random form the OneDegSq.fits catsim catalog.
The size , location and ellipticity of the second galaxy can be manually
entered by the user.
The total field is assumed to be centered at the central galaxy.
"""


def get_central_galaxy():
    """Return one randomly picked galaxy fom the catsim catalog"""
    parentdir = os.path.abspath("..")
    fname = os.path.join(parentdir, 'data', 'wldeb_data',
                         'OneDegSq.fits')
    cat = Table.read(fname, format='fits')
    np.random.seed(0)
    cond1 = (cat['fluxnorm_bulge'] != 0) & (cat['fluxnorm_disk'] != 0) \
        & (cat['fluxnorm_agn'] == 0)
    cond2 = (cat['i_ab'] < 24) & (cat['a_d'] <= 1.2) & (cat['a_d'] > 0.4)
    q, = np.where(cond1 & cond2)
    select = q[np.random.randint(0, len(q), size=1)]
    return cat[select]


def get_a_b(e, hlr):
    """Returns semimajor/minor axis from ellipticity and HLR"""
    q = (1 - e) / (1. + e)
    b = np.sqrt(hlr)
    a = b / q
    return a, b


def get_hlr(a, b):
    """Returns HLR from semimajor and minor axis"""
    hlr = (a * b)**0.5
    return hlr


def get_second_galaxy(Args, cat):
    """Assigns center and bulge, dic sizes and flux of the second galaxy
    tp the input catalog
    """
    cat['ra'][1] = cat[0]['ra'] + Args.x0 * 0.2
    cat['dec'][1] = cat[0]['dec'] + Args.y0 * 0.2
    cat['fluxnorm_bulge'][1] = cat['fluxnorm_bulge'][1] * Args.bflux_frac
    cat['fluxnorm_disk'][1] = cat['fluxnorm_disk'][1] * Args.dflux_frac
    # From ellipticity and HLR compute a, b
    hlr = get_hlr(cat['a_b'], cat['b_b'])
    a, b = get_a_b(Args.b_e, hlr * Args.bhlr_frac)
    cat['a_b'][1] = a
    cat['b_b'][1] = b
    hlr = get_hlr(cat['a_b'], cat['b_b'])
    a, b = get_a_b(Args.d_e, hlr * Args.dhlr_frac)
    cat['a_d'][1] = cat['a_d'][0] * Args.dhlr_frac**0.5
    cat['b_d'][1] = cat['b_d'][0] * Args.dhlr_frac**0.5
    cat['pa_disk'][1] = Args.p_angle
    cat['pa_bulge'][1] = Args.p_angle
    cat['galtileid'][1] = 1


def main(Args):
    # Get a catlog entry from One square degree catalog
    catalog = Table()
    cent_gal = get_central_galaxy()
    catalog = vstack([catalog, cent_gal])
    # This is the central galaxy. Assign its center to 0,0
    catalog[0]['ra'] = 0
    catalog[0]['dec'] = 0
    catalog['galtileid'][0] = 0
    # Add other galaxy
    catalog = vstack([catalog, cent_gal])
    get_second_galaxy(Args, catalog)
    parentdir = os.path.abspath("..")
    fname = os.path.join(parentdir, 'data', 'wldeb_data',
                         'gal_pair_catalog.fits')
    catalog.write(fname, format='fits', overwrite=True)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--filter_band', default="i",
                        help="LSST imaging band to simulate in[Default:i]")
    parser.add_argument('--bflux_frac', default=1,
                        help="Flux of second galaxy bulge as a fraction of \
                        central galaxy bulge flux [Default:1]")
    parser.add_argument('--dflux_frac', default=1,
                        help="Flux of second galaxy disk as a fraction of \
                        central galaxy disk flux [Default:1]")
    parser.add_argument('--bhlr_frac', default=1,
                        help="HLR of second galaxy bulge as a fraction of \
                        central galaxy bulge HLR [Default:1]")
    parser.add_argument('--dhlr_frac', default=1,
                        help="HLR of second galaxy disk as a fraction of \
                        central galaxy disk HLR [Default:1]")
    parser.add_argument('--x0', default=20,
                        help="x coordinate of center of second galaxy in pixels. \
                        Center of central galaxy is (0,0).[Default:20]")
    parser.add_argument('--y0', default=20,
                        help="y coordinate of center of second galaxy in pixels. \
                        Center of central galaxy is (0,0). [Default:20]")
    parser.add_argument('--p_angle', default=0,
                        help="Position of center of second galaxy in degrees \
                        [Default:20]")
    parser.add_argument('--b_e', default=0,
                        help="Ellipticity (e) second galaxy bulge \
                        [Default:0.2]")
    parser.add_argument('--d_e', default=0,
                        help="Ellipticity (e) second galaxy disk \
                        [Default:0.2]")
    args = parser.parse_args()
    main(args)
