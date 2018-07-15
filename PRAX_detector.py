#!/Users/cbacigalupo/canopy/bin/ python


import argparse
import logging

import socket
import sys
import traceback
import signal

class detectorClient:
    '''
    The detector client class wraps all the functionality that is needed to operate the IDL interface to the detector. 

    It provides a python interface to send the operational commands to the detector
    '''
    
    def __init__(self, host):
        self.host = host
        self.port = 5000
        self.size = 1024
        self.detectorSocket = None
        self.running = 1
        self.expTime = 1.0  # in ms
        self.loopcnt = 1
        self.frameTime =  1.47528
        self.promptSend = True
        
        self.cmdDict = {
                "ping":"PING\r\n",
                "telemetry":"GETTELEMETRY\r\n",
                "getconfig":"GETCONFIG\r\n",
                "setfsmode":"SETFSMODE(1)\r\n",
                "setramode":"SETFSMODE(0)\r\n",
                "setoutputrsf":"SETOUTPUTRESETFRAME(1)\r\n",
                "setoutputnorsf":"SETOUTPUTRESETFRAME(0)\r\n",
                "acqramp":"ACQUIRERAMP\r\n",
                # parameterized
                "setramppar":"SETRAMPPARAM",
                "setfspar":"SETFSPARAM"
                }
                
    def set_timeout(self, timeout):    
        self.detectorSocket.settimeout(timeout)
 
    def open_socket(self):
        try:
            self.detectorSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            print 'opening host:', self.host, 'port:', self.port
            self.detectorSocket.connect((self.host, self.port))
            self.detectorSocket.settimeout(30)  # .5 min
            
        except socket.error, (value, message):
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            self.terminateSocket()
            print "Could not open socket: " + message

    def terminateSocket(self):
        if self.detectorSocket:
            try:
                self.detectorSocket.close()
                self.detectorSocket.shutdown(socket.SHUT_RDWR)
                print 'Socket terminated'
            except:
                pass

    def KeyboardSignalHandler(self, signum, frame):
        print 'Keyboard Signal handler called', signum
        self.running = 0
        self.detectorSocket.close()
        self.detectorSocket.shutdown(socket.SHUT_RDWR)
        sys.exit(0)


    def SendAndReply(self, line):
        
        try:
            self.detectorSocket.send(line)
            print "line sent:", line.strip()
            
            data = self.detectorSocket.recv(self.size)
            print "line received (", len(data), "):", data.rstrip()
            
        except socket.timeout:
            print "SendAndReply: Could not write/read to/from socket. Timed out."
            data = "SendAndReply: Error"
            
        return data



    def detectorInit(self):
        
        # Open socket
        self.open_socket()
        
        # Set the signal handler 
        signal.signal(signal.SIGINT, self.KeyboardSignalHandler)
        
        print
        print 'detectorInit'
        print '------------'        

        try:
            # handle standard input
            print 'Sending ping'
            self.SendAndReply(self.cmdDict["ping"])
            print
            
        except Exception, e:
            print 'Something went wrong.', str(e)
    
    
    def getConfig(self):        
        print 'Getting config'
        self.SendAndReply(self.cmdDict["getconfig"])
        print
            
        
    def setFS(self):
        print 'Sending FS mode'
        self.SendAndReply(self.cmdDict["setfsmode"])
        print ''
        
    def setUTRamp(self):
        print 'Sending UTRamp mode'
        self.SendAndReply(self.cmdDict["setramode"])
        print ''
        
    def setFSExposure(self, nreset, nread, ngroup, exptime, nramp):
        
        print "Setting FS params..."
        self.SendAndReply("%s(%s,%s,%s,%s,%s)\r\n"
                          % (self.cmdDict["setfspar"], nreset, nread, ngroup, exptime, nramp))
        self.set_timeout(exptime*1.5)

    def setUTExposure(self, nreset, nread, ngroup, exptime, nramp):
        print "Calculating ramp params..."
        nread = 1
        ndrops = exptime/ngroup/self.frameTime-1
        print 'Using ',
        print 'nreset', nreset, 
        print 'nread', nread, 
        print 'ngroup', ngroup, 
        print 'ndrops', ndrops, 
        print 'nramp', nramp
        
        print "Setting ramp params..."
        self.SendAndReply("%s(%s,%s,%s,%s,%s)\r\n"
                          % (self.cmdDict["setramppar"], nreset, nread, ngroup, ndrops, nramp))
        self.set_timeout(exptime*1.5)
               
        
    def runExposure(self):
        
        print "Acquiring"
        self.SendAndReply(self.cmdDict["acqramp"])        
        self.set_timeout(30)




if __name__ == '__main__':

    # arg parsing for command line version 
    parser = argparse.ArgumentParser(description='PRAXIS Detector. Python wrapper to trigger exposures.')
    parser.add_argument('-host', help='Host IP. [10.80.127.5]', type=str, default='10.80.127.5', dest='host')
    parser.add_argument('-smode', help='Sampling mode. 0/1 for up the ramp or Fowler Sampling. [1]',
                         type=int, default=1, dest='smode')
    parser.add_argument('-nreset', help='Number of resets. [2]', type=int, default=2, dest='nreset')
    parser.add_argument('-nread', help='Number of reads. [1]', type=int, default=1, dest='nread')
    parser.add_argument('-ngroup', help='Number of groups [1]', type=int, default=1, dest='ngroup')
    parser.add_argument('-exptime', help='Exposure time in seconds [1]', type=float, default=1., dest='exptime')
    parser.add_argument('-nramp', help='Number of ramps. [1]', type=int, default=1, dest='nramp')
    parser.add_argument('-ndrops', help='Number of drops. [1]', type=int, default=1, dest='ndrops')
    parser.add_argument('-configonly', help='Only configuration. [True]', type=bool, default=1, dest='configOnly')
 
    args = parser.parse_args()
    print 'Command line arguments:', args
    print 

    
    
    '''
    Detector actions start here
    '''
    
    
    # initialise class and connection
    det = detectorClient(host = args.host)
    det.detectorInit()
     
    #print current config     
    det.getConfig()
    
    if not args.configOnly:
        
        # Set sampling mode
        if args.smode==0:
            det.setUTRamp()
            det.setUTExposure(nreset = args.nreset, 
                    nread = args.nread, 
                    ngroup = args.ngroup, 
                    exptime = args.exptime, 
                    nramp = args.nramp)
    
        else:
            det.setFS()
            det.setFSExposure(nreset = args.nreset, 
                            nread = args.nread, 
                            ngroup = args.ngroup, 
                            exptime = args.exptime, 
                            nramp = args.nramp)
    
    
        # expose
        det.runExposure()

        
    # wrap up
    det.terminateSocket()
