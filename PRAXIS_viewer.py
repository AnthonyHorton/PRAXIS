#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import warnings
import os
import glob

import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import matplotlib.colors as colours
from matplotlib.patches import RegularPolygon

radius_flat = 0.55 / 2  # radius to the flat edge
radius = radius_flat * 2 / 3**0.5  # radius of the circle around
radius_short = np.sin(np.pi/6) * radius  # short radius

window = ((830, 1220), (700, 1523))

tram_coef = ((63.8168, 0.00159368, -5.14315*10**-6),
             (75.6552, 0.0018184, -5.28723*10**-6),
             (87.6753, 0.00110761, -4.25759*10**-6),
             (99.742, 0.0017248, -5.03802*10**-6),
             (111.731, 0.00182788, -5.14512*10**-6),
             (123.743, 0.00134928, -4.50364*10**-6),
             (135.788, -0.000400901,  -2.13853*10**-6),
             (np.nan, np.nan, np.nan),
             (159.696, -0.00134634,  -5.88114*10**-7),
             (171.742, -0.00092513, -1.22513*10**-6),
             (183.977, -0.00178464, -1.29449*10**-6),
             (np.nan, np.nan, np.nan),
             (267.469, -0.00266766, 1.18696*10**-6),
             (279.389, -0.00292473, 1.93398*10**-6),
             (291.537, -0.00307161, 1.99748*10**-6),
             (303.404, -0.00289977,  1.7745*10**-6),
             (315.713, -0.00382432,  2.36483*10**-6),
             (327.903, -0.00511395,  3.98343*10**-6),
             (339.219, -0.00243858, 1.55742*10**-6))


def process_data(main_data, tramlines, tramlines_bg, subtract=None, throughput_file=None):
    """
    Takes image data and tramline parameters and returns the flux for each of the 19 fibres

    Args:
        filename (str) : The name of the input file
        tramlines (list of (np.array, np.array) tuples): For each fibre a tuple or arrays,
            containing the y and x coordinates of the pixels to include in the extraction.
        subtract (str, optional): filename of a sky background image to subtract from the data

    Returns:
        main_data (np.array): processed image data
        fluxes (list of floats): Summed flux for each fibre
    """
    spectra = []

    if subtract:
        subtract = check_path(subtract)
        try:
            with pf.open(subtract) as hdulist:
                sub_data = hdulist[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
        except Exception as err:
            warnings.warn('Could not open subtract file {}.'.format(args.subtract))
            raise err
        print('Subtracting {} from image\n'.format(subtract))
        main_data = main_data - sub_data
        print('Extracting spectra\n')
        for tramline in tramlines:
            spectra.append(extract_spectrum(main_data, tramline))
    elif tramlines_bg:
        print('Extracting spectra with background subtraction\n')
        for tramline, tramline_bg in zip(tramlines, tramlines_bg):
            spectra.append(extract_spectrum(main_data, tramline, tramline_bg))
    else:
        print('Extracting spectra\n')
        for tramline in tramlines:
            spectra.append(extract_spectrum(main_data, tramline))

    science_spectrum = sum(spectra[0:7])

    print('Summing flux for each fibre\n')
    fluxes = [spectrum.sum() for spectrum in spectra]

    if throughput_file:
        fluxes = throughput_correction(fluxes, throughput_file)

    for i, flux in enumerate(fluxes):
        print("{}: {}".format(i + 1, flux))
    print()

    return main_data, fluxes, science_spectrum


def make_hex_array(fluxes):
    """
    Assembles an array containing the hexagonal IFU lenslet centre coordinates together with the
    corresponding normalised fluxes.

    Args:
        fluxes (list of floats): fluxes from each fibre

    Returns:
        hex_array (np.array): hex_array[idx] = [centreX, centreY, flux]
    """
    # fibre number = idx+1
    hex_array = np.ones((19, 3)) * np.nan

    hex_array[0] = [0, 0, 0]

    hex_array[1] = [-2 * radius_flat, 0, 0]
    hex_array[2] = [-radius_flat, -(radius + radius_short), 0]
    hex_array[3] = [radius_flat, -(radius + radius_short), 0]
    hex_array[4] = [2 * radius_flat, 0, 0]
    hex_array[5] = [radius_flat, radius + radius_short, 0]
    hex_array[6] = [-radius_flat, radius + radius_short, 0]

    hex_array[7] = [-3 * radius_flat, (radius + radius_short), 0]
    hex_array[8] = [-4 * radius_flat, 0, 0]
    hex_array[9] = [-3 * radius_flat, -(radius + radius_short), 0]
    hex_array[10] = [-2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[11] = [0, -2 * (radius + radius_short), 0]
    hex_array[12] = [2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[13] = [3 * radius_flat, -(radius + radius_short), 0]
    hex_array[14] = [4 * radius_flat, 0, 0]
    hex_array[15] = [3 * radius_flat, (radius + radius_short), 0]
    hex_array[16] = [2 * radius_flat, 2 * (radius + radius_short), 0]
    hex_array[17] = [0, 2 * (radius + radius_short), 0]
    hex_array[18] = [-2 * radius_flat, 2 * (radius + radius_short), 0]

    hex_array[:, 2] = fluxes
    hex_array[:, 2] /= np.nanmax(fluxes)

    return hex_array


def initialise_tramlines(width, background_width=None):
    """
    Used tramline fit from tram_coef to produce a list of tuples of arrays, each tuple contains
    the y and x coordinates of the pixels to include in the extraction for a given fibreself. When
    formatted in this way each item in the list can be used directly to index the image data.

    Args:
        width (int): width, in pixels, of the extraction region

    Returns:
        tramlines (list of tuples of np.array): pixels coordinates for each fibre
    """
    xs = np.arange(0, window[1][1] - window[1][0])

    x_grid, y_grid = np.meshgrid(xs, np.arange(width))
    tramlines = []

    if background_width:
        x_grid_bg, y_grid_bg = np.meshgrid(xs, np.arange(background_width))
        tramlines_bg = []

    for (a, b, c) in tram_coef:
        if np.isnan(a):
            tramlines.append([])
            if background_width:
                tramlines_bg.append([])
            continue
        # Calculate curve
        ys = a + b * xs + c * xs**2
        # Calculate set of y shifted versions to get desired width
        ys_spectrum = ys.reshape((1, ys.shape[0]))
        ys_spectrum = ys_spectrum + np.arange(-(width - 1)/2, (width + 1)/2).reshape((width, 1))

        # Round to integers
        ys_spectrum = np.around(ys_spectrum, decimals=0).astype(np.int)

        # Reshape into (y coords, x coords) for numpy indexing
        tramline = (ys_spectrum.ravel(), x_grid.ravel())
        tramlines.append(tramline)

        if background_width:
            ys_bg = ys.reshape((1, ys.shape[0]))
            ys_bg = ys_bg + np.arange(-(background_width - 1)/2,
                                      (background_width + 1)/2).reshape((background_width, 1))
            ys_bg = np.around(ys_bg, decimals=0).astype(np.int)
            tramline_bg = (ys_bg.ravel(), x_grid_bg.ravel())
            tramlines_bg.append(tramline_bg)

    # Fibre number increases with decreasing y, but tram_coef is in order of increasing y.
    tramlines.reverse()

    if background_width:
        tramlines_bg.reverse()
        return tramlines, tramlines_bg

    return tramlines


def plot_hexagons(hex_array, nolabels, filename, spectrum):
    """
    Plots the reconstructed IFU image and spectrum from the brightest fibre.
    """
    fig = plt.figure(figsize=(12, 6), tight_layout=True)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect('equal')
    for i, (x, y, flux) in enumerate(hex_array):
        if np.isnan(flux):
            colour = 'red'
        else:
            colour = str(flux)
        hex = RegularPolygon((x, y),
                             numVertices=6,
                             radius=radius,
                             orientation=0,
                             facecolor=colour,
                             lw=0)
        ax1.add_patch(hex)

        # Also add a text label
        if not nolabels:
            ax1.text(x, y, i + 1, ha='center', va='center', size=10, color='green')

    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.grid(True)
    ax1.set_axisbelow(True)
    title = "{} IFU reconstruction (North up, East left)".format(filename.split('/')[-3])
    ax1.set_title(title)
    ax1.set_xlabel('Arcsecs')
    ax1.set_ylabel('Arcsecs')

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(spectrum)
    ax2.set_title('Science fibre spectrum')
    ax2.set_xlim(0, len(spectrum - 1))
    ax2.set_ylim(0, 1.05 * spectrum.max())

    plt.show()


def plot_tramlines(tramlines, image_data, tramlines_bg=None):
    """
    Displays image data with the tramline extraction regions using the viridis colour map, and
    the remainder in grey.
    """
    mask = np.ones_like(image_data)
    for tramline in tramlines:
        mask[tramline] = 0.0

    if tramlines_bg:
        mask_bg = np.ones_like(image_data)
        for tramline_bg in tramlines_bg:
            mask_bg[tramline_bg] = 0.0
        plt.imshow(np.ma.array(image_data, mask=mask_bg),
                   cmap='gray_r',
                   norm=colours.PowerNorm(gamma=0.5),
                   origin='lower')
    else:
        plt.imshow(image_data,
                   cmap='gray_r',
                   norm=colours.PowerNorm(gamma=0.5),
                   origin='lower')

    plt.imshow(np.ma.array(image_data, mask=mask),
               cmap='viridis_r',
               norm=colours.PowerNorm(gamma=0.5),
               origin='lower')
    plt.colorbar()
    plt.show()


def extract_spectrum(image_data, tramline, tramline_bg=None):
    """
    Crude, uncalibrated spectral extraction (just masks image & collapses in y direction)

    Args:
        image_data (np.array): 2D data
        tramline (tuples of np.array): y, x coordinates of the spectrum extraction region

    Returns:
        spectrum (np.array): 1D spectrum
    """
    mask = np.ones_like(image_data)
    mask[tramline] = 0.0
    spectrum = np.ma.array(image_data, mask=mask).mean(axis=0)

    if tramline_bg:
        mask_bg = np.ones_like(image_data)
        mask_bg[tramline_bg] = 0.0
        mask_bg[tramline] = 1.0
        spectrum = spectrum - np.ma.array(image_data, mask=mask_bg).mean(axis=0)

    return spectrum


def find_latest():
    prefix = os.getenv('PRAXIS')
    if not prefix:
        warnings.warn("$PRAXIS environment variable not set!")
    dirs = glob.glob((os.path.join(prefix, '20*')))
    if not dirs:
        raise FileNotFoundError("No PRAXIS data found in {}".format(prefix))
    dirs.sort()
    latest_dir = dirs[-1]
    filename = os.path.join(latest_dir, 'Result/CDSResult.fits')
    if not os.path.exists(filename):
        raise FileNotFoundError("{} contains no PRAXIS data".format(latest_dir))
    return filename


def check_path(filename):
    """
    Given a (partial) filename attempts to resolve it to a valid path to a PRAXIS data file.

    Args:
        filename (str): path, optionally a partial path, to a PRAXIS data file.

    Returns:
        filename (str): path to the PRAXIS data file.

    Notes:
        If the input filename does not end it '.fits' then 'Result/CDSResult.fits' be added. If
        the filename is still not a path to an existing file check_path will attempt to prefix
        filename with the value of the $PRAXIS environment variable.

    Raise:
        FileNotFoundError: raises if none of the above results in a path to an existing file.
    """
    _, extension = os.path.splitext(filename)
    if extension != ".fits":
        filename = os.path.join(filename, 'Result/CDSResult.fits')
    if not os.path.exists(os.path.expandvars(filename)):
        # Not a direct path to the file. Look for $PRAXIS environment variable and try
        # prefixing it with that.
        prefixed_filename = os.path.expandvars(os.path.join('$PRAXIS', filename))
        if not os.path.exists(prefixed_filename):
            # That didn't work either.
            raise FileNotFoundError("Couldn't find input file {}".format(filename))
        else:
            filename = prefixed_filename
    return filename


def throughput_correction(fluxes, throughput_file):
    """
    Reads fibre relative throughputs from a file and corrects the fibre fluxes by
    divided by the relative throughputs.

    Args:
        fluxes (np.array): initial fibre fluxes
        throughput_file (str): path to the file containing relative fibre throughput, in numpy
            text format

    Returns:
        fluxes (np.array): corrected fibre fluxes
    """
    try:
        throughputs = np.loadtxt(throughput_file)
    except Exception as err:
        warnings.warn("Error opening throughput file {}".format(throughput_file))
        raise err
    print("Correcting fibre fluxes with relative throughputs from {}\n".format(throughput_file))
    fluxes = fluxes / throughputs
    return fluxes


if __name__ == '__main__':
    # arg parsing for command line version
    description = 'PRAXIS Viewer. Visualisation tool for PRAXIS data sets.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--filename', help='Input file.', type=str)
    parser.add_argument('--width', help='Width for extraction', type=int, default=5)
    parser.add_argument('--bgwidth', help='Width background extraction', type=int, default=10)
    parser.add_argument('--tram', help='Plot of tramlines', action='store_true')
    parser.add_argument('--nolabels', help="Don't show fibre labels on plot", action='store_true')
    parser.add_argument('--subtract', help='Image for sky subtracting', type=str)
    parser.add_argument('--throughput', help='File to use for relative throughput correction',
                        type=str, default='throughput_dome_flat_20180729.txt')
    args = parser.parse_args()
    print(" *** PRAXIS Viewer *** \n")
    print("Arguments: {}\n".format(vars(args)))

    if not args.filename:
        print("Looking for latest image file...")
        filename = find_latest()
        print("Found {}\n".format(filename))
    else:
        filename = check_path(args.filename)

    try:
        with pf.open(filename) as hdulist:
            main_data = hdulist[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
    except Exception as err:
        warnings.warn('Could not open input file {}'.format(filename))
        raise err
    print('Read data from {}\n'.format(filename))

    if args.bgwidth:
        tramlines, tramlines_bg = initialise_tramlines(args.width, args.bgwidth)
    else:
        tramlines = initialise_tramlines(args.width)
        tramlines_bg = None

    if args.tram:
        plot_tramlines(tramlines, main_data, tramlines_bg)
    main_data, fluxes, science_spectrum = process_data(main_data,
                                                       tramlines=tramlines,
                                                       tramlines_bg=tramlines_bg,
                                                       subtract=args.subtract,
                                                       throughput_file=args.throughput)
    hex_array = make_hex_array(fluxes)
    if np.sum(hex_array):
        plot_hexagons(hex_array, args.nolabels, filename, science_spectrum)
