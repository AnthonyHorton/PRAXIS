#!/usr/bin/env python
import argparse
# import logging
import numpy as np
import astropy.io.fits as pf
from matplotlib import pyplot as plt
from matplotlib.patches import RegularPolygon


rf = 0.055/2  # radius to the flat edge
r = rf * 2 / 3**0.5  # radius of the circle around
rs = np.sin(np.pi/6) * r  # short radius


def parseFile(filename, tramLines):
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

    hexArray = initialiseHexArray()

    # Check if file exists
    try:
        mainFile = pf.open(filename)
    except:
        print('Could not open main file.')
        return False

    mainData = mainFile[0].data

    if args.divide:
        try:
            divFile = pf.open(args.divide)
        except:
            print('Could not open divide file.')
            return False

        print('Dividing by', args.divide)
        divData = divFile[0].data
        mainData /= divData

    if args.subtract:
        try:
            subFile = pf.open(args.subtract)
        except:
            print('Could not open subtract file.')
            return False

        subData = subFile[0].data
        if args.divide:
            subData /= divData

        print('Subtracting by', args.subtract)
        mainData -= subData

    print('Sum(flux) for each fibre:')
    for i in range(19):
        thisTramCoords = np.where(tramLines == i)
        thisFlux = np.sum(mainData[tramLines[thisTramCoords]])
        print(i+1, thisFlux)
        hexArray[i, 2] = thisFlux

    hexArray[:, 2] -= np.min(hexArray[:, 2])
    hexArray[:, 2] /= np.max(hexArray[:, 2])

    return hexArray


def initialiseHexArray():

    # hexArray[idx] = [centreX,centreY,flux]
    # fibre number = idx+1
    hexArray = np.ones((19, 3)) * np.nan

    hexArray[0] = [0, 0, 0]
    hexArray[1] = [rf, r+rs, 0]
    hexArray[2] = [2*rf, 0, 0]
    hexArray[3] = [rf, -(r+rs), 0]
    hexArray[4] = [-rf, -(r+rs), 0]
    hexArray[5] = [-2*rf, 0, 0]
    hexArray[6] = [-rf, r+rs, 0]
    hexArray[7] = [0, 2*(r+rs), 0]
    hexArray[8] = [2*rf, 2*(r+rs), 0]
    hexArray[9] = [3*rf, (r+rs), 0]
    hexArray[10] = [4*rf, 0, 0]
    hexArray[11] = [3*rf, -(r+rs), 0]
    hexArray[12] = [2*rf, -2*(r+rs), 0]
    hexArray[13] = [0, -2*(r+rs), 0]
    hexArray[14] = [-2*rf, -2*(r+rs), 0]
    hexArray[15] = [-3*rf, -(r+rs), 0]
    hexArray[16] = [-4*rf, 0, 0]
    hexArray[17] = [-3*rf, (r+rs), 0]
    hexArray[18] = [-2*rf, 2*(r+rs), 0]

    return hexArray


'''
Actually these will only make sense if you first select a subwindow:   x :700 - 1523 and y: 830 - 1218

or add 830 to a, and 700 to x.
y=a+ b x+ c xx
'''


def initialiseTramLines():
    a = []
    a.append([63.8168, 0.00159368, -5.14315*10**-6])
    a.append([75.6552, 0.0018184, -5.28723*10**-6])
    a.append([87.6753, 0.00110761, -4.25759*10**-6])
    a.append([99.742, 0.0017248, -5.03802*10**-6])
    a.append([111.731, 0.00182788, -5.14512*10**-6])
    a.append([123.743, 0.00134928, -4.50364*10**-6])
    a.append([135.788, -0.000400901,  -2.13853*10**-6])
    a.append([0., 0., 0.])
    a.append([159.696, -0.00134634,  -5.88114*10**-7])
    a.append([171.742, -0.00092513, -1.22513*10**-6])
    a.append([183.977, -0.00178464, -1.29449*10**-6])
    a.append([0., 0., 0.])
    a.append([267.469, -0.00266766, 1.18696*10**-6])
    a.append([279.389, -0.00292473, 1.93398*10**-6])
    a.append([291.537, -0.00307161, 1.99748*10**-6])
    a.append([303.404, -0.00289977,  1.7745*10**-6])
    a.append([315.713, -0.00382432,  2.36483*10**-6])
    a.append([327.903, -0.00511395,  3.98343*10**-6])
    a.append([339.219, -0.00243858, 1.55742*10**-6])
    tramCoef = np.array(a)

    tramLines = []
    #     tramCoef = np.load('tram_coef.npy')

    for i in range(19):
        a = tramCoef[i, 0]
        b = tramCoef[i, 1]
        c = tramCoef[i, 2]

        for k in range(5): #4 px with
            for x in range(700, 1523):
                y = (a+828+k + b*x + c * x**2)
                tramLines.append([i, x, y])

    tramLines = np.array(tramLines).astype(int)

    return tramLines


def plotHexagons(hexArray, args):
    '''
    Plots the fibre bundle in the spatial distribuition from an image.


    '''

    fig, ax = plt.subplots()

    ax.set_aspect('equal')
    i = 1
    for x, y, flux in zip(hexArray[:, 0], hexArray[:, 1], hexArray[:, 2]):
        thisHex = RegularPolygon((x, y),
                                 numVertices=6,
                                 radius=r,
                                 orientation=0,
                                 facecolor=str(flux),
                                 lw=0)
        ax.add_patch(thisHex)

        # Also add a text label
        if args.labels:
            ax.text(x, y, i, ha='center', va='center', size=10)
        i += 1

    plt.xlim(-0.2, 0.2)
    plt.ylim(-0.2, 0.2)
    plt.grid(True)
    ax.set_axisbelow(True)
    plt.title(args.fileName)
    plt.xlabel('Arcsecs')
    plt.ylabel('Arcsecs')
    plt.show()


def plotTramlines(filename, tramLines):

    # Check if file exists
    try:
        thisFile = pf.open(filename)
    except:
        print('Could not open file.')
        return False

    thisData = thisFile[0].data
    plt.imshow(np.log(thisData))
    plt.plot(tramLines[:, 1], tramLines[:, 2], 'r.')
    plt.show()


if __name__ == '__main__':

    # arg parsing for command line version
    parser = argparse.ArgumentParser(description='PRAXIS Viewer. Visualisation tool for PRAXIS data sets.')
    parser.add_argument('-f', help='Input file.', type=str, default='PRAXIS_sample.fits', dest='fileName')
    parser.add_argument('-tram', help='Plot of tramlines', type=bool, default=False, dest='tram')
    parser.add_argument('-labels', help='Show fibre labels on plot', type=bool, default=False, dest='labels')
    parser.add_argument('-divide', help='Image for flat dividing', type=str, dest='divide')
    parser.add_argument('-subtract', help='Image for sky subtracting', type=str, dest='subtract')
    args = parser.parse_args()

    print(args)
    print()

    tramLines = initialiseTramLines()

    if args.tram:
        plotTramlines(args.fileName, tramLines)

    hexArray = parseFile(args.fileName, tramLines)

    print()
    print('HexArray (x,y, Normalised flux):')
    print(hexArray)
    if np.sum(hexArray):
        plotHexagons(hexArray, args)
