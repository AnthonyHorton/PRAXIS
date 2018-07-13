
import argparse
import logging
import numpy as np




if __name__ == '__main__':

    #arg parsing for command line version 
    parser = argparse.ArgumentParser(description='PRAXIS Viewer. Visualisation tool for PRAXIS data sets.')
    parser.add_argument('-v', help='Verbosity level (0-None, 5-Max).')  
    parser.add_argument('-o', help='Routed tile output file name (S3)', type=str, metavar='RTileFileNameS3.json')
    parser.add_argument('-f', help='Allocation target file (S2)', type=str, metavar='XYTileFileNameS2.json', default='s2_example.json')
#     parser.add_argument('-fcsv', help='filename of the output csv file', default='out.csv')
#     parser.add_argument('-fjson', help='filename of the output json file', default='out.json')
    args = parser.parse_args()
    print args


def parseFile():
    '''
    Takes an input file and tramline parameters and returns the flux for each of the 19 fibres
    
    Args:
        filename (str) : The name of the input file
        tramLines (np.ndarray) : 3D - [[x,y],..],...] list of pixels to sum over for each fibre
        
    Returns:
        hexArray(np.ndarray) = 2D - [[centreX,centreY,flatRad,flux],...] 
        
    Note:
        Creating random data at the moment until tramlines are calculated
  
    '''
    r = 0.055  # flat radius
    hexArray = np.ones((19,4))
    
    hexArray[0] = [0,0,r,0]
    hexArray[1] = [0,0,0.055,0]
    hexArray[2] = [0,0,0.055,0]
    hexArray[3] = [0,0,0.055,0]
    hexArray[4] = [0,0,0.055,0]
    hexArray[5] = [0,0,0.055,0]
    hexArray[6] = [0,0,0.055,0]
    hexArray[7] = [0,0,0.055,0]
    hexArray[8] = [0,0,0.055,0]
    hexArray[9] = [0,0,0.055,0]
    hexArray[10] = [0,0,0.055,0]
    hexArray[11] = [0,0,0.055,0]
    hexArray[12] = [0,0,0.055,0]
    hexArray[13] = [0,0,0.055,0]
    hexArray[14] = [0,0,0.055,0]
    hexArray[15] = [0,0,0.055,0]
    hexArray[16] = [0,0,0.055,0]
    hexArray[17] = [0,0,0.055,0]
    hexArray[18] = [0,0,0.055,0]
