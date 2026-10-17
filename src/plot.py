description = """
Small command line application for using plot_viewer.py to plot/visualize data from generated .FITS files and to run 
data file generation scripts.
"""

import sys
import argparse
from daomop import storage
from plotting import plot_viewer, st_deviation, mag_data, latidutes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--std-histogram", "-sh", dest="std_hist", default=False, action="store_true",
                        help="PLOT: Create magnitude standard deviation histograms.")

    parser.add_argument("--object-count", "-oc", dest="lat_count", default=False, action="store_true",
                        help="PLOT: Create object count per latitude histograms.")

    parser.add_argument("--sky-coverage", "-sc", dest="sky_cov", default=False, action="store_true",
                        help="PLOT: Create sky coverage data scatter plots.")

    parser.add_argument("--st-deviation", "-std", dest="std", default=False, action="store_true",
                        help="Run st_deviation_plot.py to generate standard deviation data files.")

    parser.add_argument("--latitudes", "-lat", dest="lat", default=False, action="store_true",
                        help="Run latitudes.py to generate object count per latitude data files.")

    parser.add_argument("--magnitudes", "-mag", dest="mag", default=False, action="store_true",
                        help="Run mag_data.py to generate coverage data files.")

    args = parser.parse_args()

    if len(sys.argv[1:]) > 1:
        parser.error("Too many arguments: expected exactly one. Use -h for help.")

    elif args.std_hist or args.lat_count or args.sky_cov:
        plot_viewer.main(args)
    elif args.std:
        st_deviation.main()
    elif args.mag:
        mag_data.main()
    elif args.lat:
        latidutes.main()
    else:
        parser.error('No action requested, add exactly one argument (use -h to see shorthand notation)')

__doc__ = description
