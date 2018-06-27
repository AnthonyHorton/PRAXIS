
import argparse
import logging





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


