"""Script creates pair of galaxy with input parametrs"""
import argparse


def main():
    """makes a catalog that can be input to WLdeblending package to create
    image with 2 galaxy blend.
    """


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--filter_band', default="i",
                        help="LSST imaging band to simulate in[Default:i]")
    parser.add_argument('--bflux_frac', default=1,
                        help="Fraction of central galaxy bulge flux assigned \
                        to second galaxy [Default:1]")
    parser.add_argument('--dflux_frac', default=1,
                        help="Fraction of central galaxy bulge flux assigned \
                        to second galaxy [Default:1]")
    args = parser.parse_args()
    main()
