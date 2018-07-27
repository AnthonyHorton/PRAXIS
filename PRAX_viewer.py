#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import warnings
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import matplotlib.colors as colours
from matplotlib.patches import RegularPolygon


radius_flat = 0.055 / 2  # radius to the flat edge
radius = radius_flat * 2 / 3**0.5  # radius of the circle around
radius_short = np.sin(np.pi/6) * radius  # short radius

window = ((828, 1218), (700, 1523))

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


def parse_file(filename, tramlines):
    '''
    Takes an input file and tramline parameters and returns the flux for each of the 19 fibres

    Args:
        filename (str) : The name of the input file
        tramLines (np.ndarray) : 3D - [[x,y],..],...] list of pixels to sum over for each fibre

    Returns:
        hexArray(np.ndarray) = 2D - [[centreX,centreY,flux],...]

    Note:
        Creating random data at the moment until tramlines are calculated

    '''

    hex_array = initialise_hex_array()

    # Check if file exists
    try:
        main_file = pf.open(filename)
    except Exception as err:
        warnings.warn('Could not open main file {}'.format(filename))
        raise err

    main_data = main_file[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
    if args.divide:
        try:
            div_file = pf.open(args.divide)
        except Exception as err:
            warnings.warn('Could not open divide file {}.'.format(args.divide))
            raise err

        print('Dividing {} by {}'.format(fileame, arg.divide))
        div_data = div_file[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
        main_data /= div_data

    if args.subtract:
        try:
            sub_file = pf.open(args.subtract)
        except Exception as err:
            warnings.warn('Could not open subtract file {}.'.format(args.subtract))
            raise err

        sub_data = sub_file[0].data[window[0][0]:window[0][1], window[1][0]:window[1][1]]
        if args.divide:
            sub_data /= div_data

        print('Subtracting {} from {}'.format(args.subtract, filename))
        main_data -= sub_data

    print('Sum flux for each fibre:')
    for i, tramline in enumerate(tramlines):
        flux = np.sum(main_data[tramline])
        print(i+1, flux)
        if flux == 0:
            hex_array[i, 2] = np.nan
        else:
            hex_array[i, 2] = flux

    hex_array[:, 2] -= np.nanmin(hex_array[:, 2])
    hex_array[:, 2] /= np.nanmax(hex_array[:, 2])

    return hex_array


def initialise_hex_array():

    # hexArray[idx] = [centreX,centreY,flux]
    # fibre number = idx+1
    hex_array = np.ones((19, 3)) * np.nan

    hex_array[0] = [0, 0, 0]
    hex_array[1] = [radius_flat, radius + radius_short, 0]
    hex_array[2] = [2 * radius_flat, 0, 0]
    hex_array[3] = [radius_flat, -(radius + radius_short), 0]
    hex_array[4] = [-radius_flat, -(radius + radius_short), 0]
    hex_array[5] = [-2 * radius_flat, 0, 0]
    hex_array[6] = [-radius_flat, radius + radius_short, 0]
    hex_array[7] = [0, 2 * (radius + radius_short), 0]
    hex_array[8] = [2 * radius_flat, 2 * (radius + radius_short), 0]
    hex_array[9] = [3 * radius_flat, (radius + radius_short), 0]
    hex_array[10] = [4 * radius_flat, 0, 0]
    hex_array[11] = [3 * radius_flat, -(radius + radius_short), 0]
    hex_array[12] = [2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[13] = [0, -2 * (radius + radius_short), 0]
    hex_array[14] = [-2 * radius_flat, -2 * (radius + radius_short), 0]
    hex_array[15] = [-3 * radius_flat, -(radius + radius_short), 0]
    hex_array[16] = [-4 * radius_flat, 0, 0]
    hex_array[17] = [-3 * radius_flat, (radius + radius_short), 0]
    hex_array[18] = [-2 * radius_flat, 2 * (radius + radius_short), 0]

    return hex_array


'''
Actually these will only make sense if you first select a subwindow:   x :700 - 1523 and y: 830 - 1218

or add 830 to a, and 700 to x.
y=a+ b x+ c xx
'''

def initialise_tramlines(width=4):

    xs = np.arange(0, window[1][1] - window[1][0])
    x_grid, y_grid = np.meshgrid(xs, np.arange(width))

    tramlines = []

    for (a, b, c) in tram_coef:
        if np.isnan(a):
            tramlines.append([])
            continue
        # Calculate curve
        ys = a + b * xs + c * xs**2
        # Calculate set of y shifted versions to get desired width
        ys = ys.reshape((1, ys.shape[0]))
        ys = ys + np.arange(-(width - 1)/2, (width + 1)/2).reshape((width, 1))
        # Round to integers
        ys = np.around(ys, decimals=0).astype(np.int)
        # Reshape into (y coords, x coords) for numpy indexing
        tramline = (ys.ravel(), x_grid.ravel())
        tramlines.append(tramline)

    return tramlines


def plot_hexagons(hex_array, args):
    '''
    Plots the fibre bundle in the spatial distribuition from an image.

    '''

    fig, ax = plt.subplots()

    ax.set_aspect('equal')
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
        ax.add_patch(hex)

        # Also add a text label
        if args.labels:
            ax.text(x, y, i + 1, ha='center', va='center', size=10, color='green')

    plt.xlim(-0.2, 0.2)
    plt.ylim(-0.2, 0.2)
    plt.grid(True)
    ax.set_axisbelow(True)
    plt.title(args.filename)
    plt.xlabel('Arcsecs')
    plt.ylabel('Arcsecs')
    plt.show()


def plot_tramlines(filename, tramlines):

    # Check if file exists
    try:
        image_data = pf.getdata(filename)[window[0][0]:window[0][1], window[1][0]:window[1][1]]
    except Exception as err:
        warnings.warn('Could not open file {}'.format(filename))
        raise err

    mask = np.ones_like(image_data)
    for tramline in tramlines:
        mask[tramline] = 0.0

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


if __name__ == '__main__':

    # arg parsing for command line version
    parser = argparse.ArgumentParser(description='PRAXIS Viewer. Visualisation tool for PRAXIS data sets.')
    parser.add_argument('filename', help='Input file.', type=str)
    parser.add_argument('--width', help='Width for extraction', type=int)
    parser.add_argument('--tram', help='Plot of tramlines', action='store_true')
    parser.add_argument('--labels', help='Show fibre labels on plot', action='store_true')
    parser.add_argument('--divide', help='Image for flat dividing', type=str)
    parser.add_argument('--subtract', help='Image for sky subtracting', type=str)
    args = parser.parse_args()

    print(args)
    print()

    tramlines = initialise_tramlines(args.width)

    if args.tram:
        plot_tramlines(args.filename, tramlines)

    hex_array = parse_file(args.filename, tramlines)

    print()
    print('Hex_array (x, y, normalised flux):')
    print(hex_array)
    if np.sum(hex_array):
        plot_hexagons(hex_array, args)
