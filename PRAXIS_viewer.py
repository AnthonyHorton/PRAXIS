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

window = ((828, 1218), (700, 1523))


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


def initialise_tramlines(tram_coef_file, width, background_width=None):
    """
    Used tramline fit from tram_coef to produce a list of tuples of arrays, each tuple contains
    the y and x coordinates of the pixels to include in the extraction for a given fibre. When
    formatted in this way each item in the list can be used directly to index the image data.

    Args:
        width (int): width, in pixels, of the spectral extraction region
        background_width (int): total width, in pixels, of the spectral extraction regions including
            the instrument background subtraction regions.

    Returns:
        tramlines (list of tuples of np.array): pixels coordinates for each fibre
        tramlines_bg (list of tuples of np.array): pixel coordinates for each fibre, including
            background regions.
    """
    try:
        tram_coef = np.loadtxt(tram_coef_file, delimiter=',')
    except Exception as err:
        warnings.warn("Couldn't read tramline coefficients file {}".format(tram_coef_file))
        raise err

    print("Initialising tramlines with coefficients from {}\n".format(tram_coef_file))

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


def extract_spectrum(image_data, tramline, tramline_bg=None):
    """
    Crude, uncalibrated spectral extraction (just masks image & collapses in y direction)

    Args:
        image_data (np.array): 2D data
        tramline (tuples of np.array): y, x coordinates of the spectrum extraction region
        tramline_bg (tuples of np.array, optional):

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

    return spectrum.filled(fill_value=np.nan)


def plot_tramlines(image_data, tramlines, tramlines_bg=None):
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
    hex_array[1] = [-radius_flat, radius + radius_short, 0]
    hex_array[2] = [-2 * radius_flat, 0, 0]
    hex_array[3] = [-radius_flat, -(radius + radius_short), 0]
    hex_array[4] = [radius_flat, -(radius + radius_short), 0]
    hex_array[5] = [2 * radius_flat, 0, 0]
    hex_array[6] = [radius_flat, radius + radius_short, 0]
    hex_array[7] = [0, 2 * (radius + radius_short), 0]
    hex_array[8] = [-2 * radius_flat, 2 * (radius + radius_short), 0]
    hex_array[9] = [-3 * radius_flat, (radius + radius_short), 0]
    hex_array[10] = [-4 * radius_flat, 0, 0]
    hex_array[11] = [-3 * radius_flat, -(radius + radius_short), 0]
    hex_array[12] = [-2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[13] = [0, -2 * (radius + radius_short), 0]
    hex_array[14] = [2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[15] = [3 * radius_flat, -(radius + radius_short), 0]
    hex_array[16] = [4 * radius_flat, 0, 0]
    hex_array[17] = [3 * radius_flat, (radius + radius_short), 0]
    hex_array[18] = [2 * radius_flat, 2 * (radius + radius_short), 0]

    hex_array[:, 2] = fluxes
    hex_array[:, 2] /= np.nanmax(fluxes)

    return hex_array


def plot_hexagons(hex_array, filename, spectrum):
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


def process_data(filename,
                 subtract,
                 tram_coeffs,
                 width,
                 background_width,
                 plot_tram,
                 throughput_file,
                 plot_hex):
    """
    Takes image data and tramline parameters and returns the flux for each of the 19 fibres

    Args:
        filename (str, optional) : the name of the input file. If not given will try to find the
            latest PRAXIS file and process that. If given a partial path (e.g. just the datetime
            part) it will attempt to autocomplete.
        subtract (str, optional): filename of a sky background image to subtract from the data. If
            not given instrument background subtraction will be used instead.
        tram_coeffs (str, optional): filename of the spectum tramline coefficients.
        width (int): Width, in pixels, of the spectral extraction regions
        background_width (int): total width, in pixels, of the spectral extraction regions including
            the instrument background subtraction regions
        plot_tram (bool): if True will display the image data with tramlines overlaid.
        throughput_file (str): filename of the fibre relative throughputs
        plot_hex (bool): if True will display the IFU image reconstruction and science spectrum

    Returns:
        main_data (np.array): processed image data
        fluxes (list of floats): Summed flux for each fibre
        science_spectrum (np.array): Summed spectrum from the 7 science fibres
    """
    # Find and open main data file
    if filename is None:
        print("Looking for latest image file...")
        filename = find_latest()
        print("Found {}\n".format(filename))
    else:
        filename = check_path(filename)
    try:
        with pf.open(filename) as hdulist:
            main_data = hdulist[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
    except Exception as err:
        warnings.warn('Could not open input file {}'.format(filename))
        raise err
    print('Read data from {}\n'.format(filename))

    # Background subtraction and spectral extraction
    spectra = []
    if not subtract:
        tramlines, tramlines_bg = initialise_tramlines(tram_coeffs, width, background_width)
        print('Extracting spectra with background subtraction\n')
        for tramline, tramline_bg in zip(tramlines, tramlines_bg):
            spectra.append(extract_spectrum(main_data, tramline, tramline_bg))
    else:
        tramlines = initialise_tramlines(tram_coeffs, width)
        tramlines_bg = None
        subtract = check_path(subtract)
        try:
            with pf.open(subtract) as hdulist:
                sub_data = hdulist[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
        except Exception as err:
            warnings.warn('Could not read subtract file {}.'.format(args.subtract))
            raise err
        print('Subtracting {} from image\n'.format(subtract))
        main_data = main_data - sub_data
        print('Extracting spectra\n')
        for tramline in tramlines:
            spectra.append(extract_spectrum(main_data, tramline))

    if plot_tram:
        plot_tramlines(tramlines, main_data, tramlines_bg)

    science_spectrum = np.nanmean(np.array(spectra[0:7]), axis=0)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Mean of empty slice')
        fluxes = [np.nanmean(spectrum) for spectrum in spectra]
    if throughput_file:
        fluxes = throughput_correction(fluxes, throughput_file)

    print("Fibre fluxes:\n")
    for i, flux in enumerate(fluxes):
        print("{}: {}".format(i + 1, flux))
    print()

    if plot_hex:
        hex_array = make_hex_array(fluxes)
        if np.sum(hex_array):
            plot_hexagons(hex_array, filename, science_spectrum)

    return main_data, fluxes, science_spectrum


if __name__ == '__main__':
    # arg parsing for command line version
    description = 'PRAXIS Viewer. Visualisation tool for PRAXIS data sets.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--filename', help='Input file.', type=str)
    parser.add_argument('--subtract', help='Sky subtraction file', type=str)
    parser.add_argument('--tramcoeffs', help='Tramline coefficients file', type=str,
                        default='traces.csv')
    parser.add_argument('--width', help='Width for spectral extraction', type=int, default=5)
    parser.add_argument('--bgwidth', help='Width for background extraction', type=int, default=10)
    parser.add_argument('--plottram', help='Plot tramlines', action='store_true')
    parser.add_argument('--throughputs', help='Fibre throughputs file', type=str,
                        default='throughput_dome_flat_20180729.txt')
    parser.add_argument('--plothex', help='Plot IFU image', action='store_true', default=True)
    args = parser.parse_args()

    print(" *** PRAXIS Viewer *** \n")
    print("Arguments: {}\n".format(vars(args)))

    main_data, fluxes, science_spectrum = process_data(filename=args.filename,
                                                       subtract=args.subtract,
                                                       tram_coeffs=args.tramcoeffs,
                                                       width=args.width,
                                                       background_width=args.bgwidth,
                                                       plot_tram=args.plottram,
                                                       throughput_file=args.throughputs,
                                                       plot_hex=args.plothex)
