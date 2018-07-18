
import argparse
# import logging
import numpy as np
import pyfits as pf
import pylab as plt
from matplotlib.patches import RegularPolygon


rf = 0.055/2 # radius to the flat edge
r = rf * 2 /3**0.5 # radius of the circle around 
rs = np.sin(np.pi/6) * r # short radius



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
    
    #Check if file exists
    try:
        mainFile = pf.open(filename)
    except: 
        print 'Could not open main file.'
        return False
    
    
    mainData = mainFile[0].data
    
    if args.divide:
        try:
            divFile = pf.open(args.divide)
        except: 
            print 'Could not open divide file.'
            return False

        print 'Dividing by', args.divide
        divData = divFile[0].data
        mainData /= divData
        
    if args.subtract:
        try:
            subFile = pf.open(args.subtract)
        except: 
            print 'Could not open subtract file.'
            return False
        
        subData = subFile[0].data
        if args.divide:
            subData /= divData
        
        print 'Subtracting by', args.subtract
        mainData -= subData
        
    
    print 'Sum(flux) for each fibre:'
    for i in range(19):
        thisTramCoords = np.where(tramLines==i)
        thisFlux = np.sum(mainData[tramLines[thisTramCoords]])
        print i+1, thisFlux
        hexArray[i,2] = thisFlux
    
    hexArray[:,2] -= np.min(hexArray[:,2])
    hexArray[:,2] /= np.max(hexArray[:,2])
    
    return hexArray

def initialiseHexArray():
    
    
    #hexArray[idx] = [centreX,centreY,flux]
    #fibre number = idx+1
    hexArray = np.ones((19,3)) *np.nan
    
    hexArray[0] = [0,0,0]
    hexArray[1] = [rf,r+rs,0]
    hexArray[2] = [2*rf,0,0]
    hexArray[3] = [rf,-(r+rs),0]
    hexArray[4] = [-rf,-(r+rs),0]
    hexArray[5] = [-2*rf,0,0]
    hexArray[6] = [-rf,r+rs,0]
    hexArray[7] = [0,2*(r+rs),0]
    hexArray[8] = [2*rf,2*(r+rs),0]
    hexArray[9] = [3*rf,(r+rs),0]
    hexArray[10] = [4*rf,0,0]
    hexArray[11] = [3*rf,-(r+rs),0]
    hexArray[12] = [2*rf,-2*(r+rs),0]
    hexArray[13] = [0,-2*(r+rs),0]
    hexArray[14] = [-2*rf,-2*(r+rs),0]
    hexArray[15] = [-3*rf,-(r+rs),0]
    hexArray[16] = [-4*rf,0,0]
    hexArray[17] = [-3*rf,(r+rs),0]
    hexArray[18] = [-2*rf,2*(r+rs),0]

    return hexArray

def initialiseTramLines():
    
    tramLines = []
    
    for i in range(19):
        for j in range(500, 1500):
            kStart = 1068+i*50
            kEnd = kStart + 4
            for k in range(kStart, kEnd):
                tramLines.append([i,j,k])

    tramLines = np.array(tramLines)

    return tramLines


def plotHexagons(hexArray, args):
    '''
    Plots the fibre bundle in the spatial distribuition from an image. 
    
    
    '''
    
    fig, ax = plt.subplots()
    
    ax.set_aspect('equal')
    i=1
    for x, y, flux in zip(hexArray[:,0], hexArray[:,1], hexArray[:,2]):
        thisHex = RegularPolygon((x, y), numVertices=6, radius=r, 
                             orientation=0, 
                             facecolor=str(flux), lw=0)
        ax.add_patch(thisHex)
    
        # Also add a text label
        if args.labels: ax.text(x, y, i, ha='center', va='center', size=10)
        i+=1
        
    plt.xlim(-0.2,0.2)
    plt.ylim(-0.2,0.2)
    plt.grid(True)
    ax.set_axisbelow(True)
    plt.title(args.fileName)
    plt.xlabel('Arcsecs')
    plt.ylabel('Arcsecs')
    plt.show()
    
    
def plotTramlines(filename, tramLines):
    
    #Check if file exists
    try:
        thisFile = pf.open(filename)
    except: 
        print 'Could not open file.'
        return False

    thisData = thisFile[0].data
    plt.imshow(thisData)
    plt.plot(tramLines[:,1],tramLines[:,2], 'r.')
    plt.show()
    
if __name__ == '__main__':

    #arg parsing for command line version 
    parser = argparse.ArgumentParser(description='PRAXIS Viewer. Visualisation tool for PRAXIS data sets.')
    parser.add_argument('-f', help='Input file.', type=str, default='PRAXIS_sample.fits', dest='fileName')
    parser.add_argument('-tram', help='Plot of tramlines', type=bool, default=False, dest='tram')
    parser.add_argument('-labels', help='Show fibre labels on plot', type=bool, default=False, dest='labels')
    parser.add_argument('-divide', help='Image for flat dividing', type=str, dest='divide')
    parser.add_argument('-subtract', help='Image for sky subtracting', type=str, dest='subtract')
    args = parser.parse_args()
    
    print args
    print

    tramLines = initialiseTramLines()
    
    if args.tram:
        plotTramlines(args.fileName, tramLines)
        
    hexArray  = parseFile(args.fileName, tramLines)
    
    print
    print 'HexArray (x,y, Normalised flux):'
    print hexArray
    if np.sum(hexArray):
        plotHexagons(hexArray, args)
        
        
    
