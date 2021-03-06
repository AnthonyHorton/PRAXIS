#!/usr/bin/env python
from __future__ import print_function, division
import argparse

from praxis import process_data

if __name__ == '__main__':
    # arg parsing for command line version
    description = 'PRAXIS Viewer. Visualisation tool for PRAXIS data sets.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--filenames', help='Input file(s).', nargs='*', type=str)
    parser.add_argument('--subtract', help='Sky subtraction file(s)', nargs='*', type=str)
    parser.add_argument('--pixelmask', help='Bad pixel mask file', type=str,
                        default='calibration_data/bad_pixels/bad_pixels_20190716.txt')
    parser.add_argument('--tramcoeffs', help='Tramline coefficients file', type=str,
                        default='calibration_data/traces/traces_20190717.csv')
    parser.add_argument('--width', help='Width for spectral extraction', type=int,
                        default=5)
    parser.add_argument('--bgwidth', help='Width for background extraction', type=int,
                        default=10)
    parser.add_argument('--plottram', help='Plot tramlines', action='store_true')
    parser.add_argument('--throughputs', help='Fibre throughputs file', type=str,
                        default='calibration_data/throughputs/throughput_dome_flat_20190717.txt')
    parser.add_argument('--plothex', help='Plot IFU image', action='store_true',
                        default=True)
    parser.add_argument('--sigmaclip', help='Threshold for sigma clipping of science spectra.',
                        type=float)
    parser.add_argument('--standard', help='Divide science spectrum buy standard spectrum',
                        action='store_true', default=False)
    parser.add_argument('--standardfile', help='Standard star spectrum to correct science spectrum',
                        type=str,
                        default='calibration_data/standards/suppressed_20190718_HIP87881.txt')
    args = parser.parse_args()

    print(" *** PRAXIS Viewer *** \n")
    print("Arguments: {}\n".format(vars(args)))

    main_data, fluxes, science_spectrum = process_data(filenames=args.filenames,
                                                       subtract=args.subtract,
                                                       pixel_mask=args.pixelmask,
                                                       tram_coeffs=args.tramcoeffs,
                                                       width=args.width,
                                                       background_width=args.bgwidth,
                                                       plot_tram=args.plottram,
                                                       throughput_file=args.throughputs,
                                                       plot_hex=args.plothex,
                                                       sigma_clip=args.sigmaclip,
                                                       standard=args.standard,
                                                       standard_file=args.standardfile)
