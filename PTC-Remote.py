# -*- coding: utf-8 -*-
"""
Created on Wed Aug 09 09:31:47 2017

@author: eloy
"""

#!/usr/bin/env python

"""
T E S T H2RG
acquire ptc images
Version 0.1 2017-08-03 TFe
"""

"""
 uses 
"""
#import subprocess
#import select as uselect
import socket
import sys
#import os
import traceback
#import threading
#import array
#import string
import signal
#import time
#import struct
#import getopt
#import pyfits # load FITS module
#import tty
#import termios
import platform
if platform.system() == "Windows":
   import msvcrt

#from numpy import *

class ptc:
    def __init__(self):
        # TFe
        # local settings are hostAzCam = IPAddr, hostImage = localhost
        #
        self.hostH2RG = '141.33.128.169'
        self.portH2RG = 5000
        self.portImage = 5001
        self.size = 1024
        self.serverH2RG = None
        self.threads = None
        self.running = 1
        self.expTime = 1.0 # in ms
        self.loopcnt = 1
        self.promptSend = True
        self.cmds = {
                "ping":"PING\r\n",
                "telemetry":"GETTELEMETRY\r\n",
                "getconfig":"GETCONFIG\r\n",
                "setfsmode":"SETFSMODE(1)\r\n",
                "setramode":"SETFSMODE(0)\r\n",
                "setoutputrsf":"SETOUTPUTRESETFRAME(1)\r\n",
                "setoutputnorsf":"SETOUTPUTRESETFRAME(0)\r\n",
                "acqramp":"ACQUIRERAMP\r\n",
                # parameterized
                "setfspar":"SETFSPARAM"
                }
                
    def set_timeout(self,timeout):    
        self.serverH2RG.settimeout(timeout) # 15 min
 
    def open_socket(self):
        try:
            self.serverH2RG = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            print "opening host: "+self.hostH2RG+" port: "+`self.portH2RG`
            self.serverH2RG.connect((self.hostH2RG,self.portH2RG))
            self.serverH2RG.settimeout(30) # .5 min
        except socket.error, (value,message):
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            if self.serverH2RG:
                self.serverH2RG.close()
            print "Could not open socket: " + message
            print "*** lineno:", traceback.tb_lineno(exceptionTraceback)
            sys.exit(1)

    def terminate(self):
        if self.serverH2RG:
            try:
                self.serverH2RG.close()
                self.serverH2RG.shutdown(socket.SHUT_RDWR)
            except:
                pass

    def KeyboardSignalHandler(self, signum, frame):
        print 'Keyboard Signal handler called', signum
        self.running = 0
        self.serverH2RG.close()
        self.serverH2RG.shutdown(socket.SHUT_RDWR)
        #self.threads.stop()
        sys.exit(0)

    def SendAndReply(self, line):
        try:
            self.serverH2RG.send(line)
            print "line sent: "+line
            data = self.serverH2RG.recv(self.size)
            print "received = "+`len(data)`+" >>> "+data.rstrip()+" <<<"
            #sys.stdout.write(data)
        except socket.timeout:
            print "SendAndReply: Could not write/read to/from socket"
            data = "SendAndReply: Error"
        return data

    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False


    def setFSParam(self,nreset,nread,ngroup,exptime,nramp):
           fspar = self.SendAndReply("%s(%s,%s,%s,%s,%s)\r\n"%
                   (self.cmds["setfspar"],nreset,nread,ngroup,exptime,nramp))
           print fspar

    def ptc_init(self):
        self.open_socket()
        # Set the signal handler 
        signal.signal(signal.SIGINT, self.KeyboardSignalHandler)

        print "ptc_init: "

        try:
           # handle standard input

           ping = self.SendAndReply(self.cmds["ping"])
           print ping
           #getcfg = self.SendAndReply(self.cmds["getcfg"])
           #print getcfg
           fsmode = self.SendAndReply(self.cmds["setfsmode"])
           print fsmode

        except:
           pass
    def run(self):
        print "setting params"
        self.setFSParam("2","32","1","2400","1")
        print "acquire..."
        self.set_timeout(2700)
        fsacq = self.SendAndReply(self.cmds["acqramp"])
        print fsacq
        self.set_timeout(30)

if __name__ == "__main__":
    p = ptc()
    p.ptc_init()
    
    for ii in range(25):
        print "Sample number " + str(ii)
        p.run()
        print "Acquired number " + str(ii)
        
    #p.run(10)
    p.terminate()
    sys.exit(0)

# number of reads is set to 2
# Always set time out larger than the longest exposure, if not error of connection
# Longest exposure*1.5 = timeout
