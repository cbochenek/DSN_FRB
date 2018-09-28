#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 00:18:37 2017

@author: Aaron B. Pearlman
         aaron.b.pearlman@caltech.edu
         Division of Physics, Mathematics, and Astronomy
         California Institute of Technology
         Jet Propulsion Laboratory

m_convertPRESTODatFileToFilterbank.py - Read in a PRESTO .dat file and the
                                        information from its corresponding .inf
                                        file. Convert the data into filterbank
                                        format.

"""

import numpy as np
import sys
import os
import struct
from subprocess import call, check_output

import filterbank_convert as fil
import filterbank_class as fil2

import matplotlib.pyplot as plt



""" Read data from the PRESTO .dat file. """
def readPrestoDatFile(prestoDatFilename):

    prestoData = np.fromfile(prestoDatFilename, dtype="float32")
    prestoSamples = np.arange(1, len(prestoData) + 1, 1)

    return prestoSamples, prestoData;



""" Structure for storing data read in from the *.inf file. """
class infDataStruct:
    def __init__(self, dataFileNameString, telescopeNameString, instrumentString, \
    objectNameString, RAString, DECString, dataObservedByString, epochTimeMJDString, \
    barycenterFlagString, nbinsTimeSeriesString, timeResolutionString, \
    dataBreakFlagString, EMBandString, beamDiameterString, dispersionMeasureString, \
    centralFreqBottomBandString, bandwidthString, numChannelsString, \
    channelBandwidthString, dataAnalyzedByString):

        self.dataFileNameString = dataFileNameString
        self.telescopeNameString = telescopeNameString
        self.instrumentString = instrumentString
        self.objectNameString = objectNameString
        self.RAString = RAString
        self.DECString = DECString
        self.dataObservedByString = dataObservedByString
        self.epochTimeMJDString = epochTimeMJDString
        self.barycenterFlagString = barycenterFlagString
        self.nbinsTimeSeriesString = nbinsTimeSeriesString
        self.timeResolutionString = timeResolutionString
        self.dataBreakFlagString = dataBreakFlagString
        self.EMBandString = EMBandString
        self.beamDiameterString = beamDiameterString
        self.dispersionMeasureString = dispersionMeasureString
        self.centralFreqBottomBandString = centralFreqBottomBandString
        self.bandwidthString = bandwidthString
        self.numChannelsString = numChannelsString
        self.channelBandwidthString = channelBandwidthString
        self.dataAnalyzedByString = dataAnalyzedByString



""" Read the information fro the PRESTO .inf file. """
def readPrestoInfFile(infFileName):

    readFile = open(infFileName, 'r')

    # Read the information from the .inf file.
    dataFileNameReadLine = readFile.readline().strip()
    telescopeNameReadLine = readFile.readline().strip()
    instrumentReadLine = readFile.readline().strip()
    objectNameReadLine = readFile.readline().strip()
    RAReadLine = readFile.readline().strip()
    DECReadLine = readFile.readline().strip()
    dataObservedByReadLine = readFile.readline().strip()
    epochTimeMJDReadLine = readFile.readline().strip()
    barycenterFlagReadLine = readFile.readline().strip()
    nbinsTimeSeriesReadLine = readFile.readline().strip()
    timeResolutionReadLine = readFile.readline().strip()
    dataBreakFlagReadLine = readFile.readline().strip()
    EMBandReadLine = readFile.readline().strip()
    beamDiameterReadLine = readFile.readline().strip()
    dispersionMeasureReadLine = readFile.readline().strip()
    centralFreqBottomBandReadLine = readFile.readline().strip()
    bandwidthReadLine = readFile.readline().strip()
    numChannelsReadLine = readFile.readline().strip()
    channelBandwidthReadLine = readFile.readline().strip()
    dataAnalyzedByReadLine = readFile.readline().strip()

    # Set the delimiter
    delimiter = "= "


    # Parse the strings from the control file.
    dataFileNameString = (" ".join(dataFileNameReadLine.split())).split(delimiter)
    dataFileNameString = dataFileNameString[1]

    telescopeNameString = (" ".join(telescopeNameReadLine.split())).split(delimiter)
    telescopeNameString = telescopeNameString[1]

    instrumentString = (" ".join(instrumentReadLine.split())).split(delimiter)
    instrumentString = instrumentString[1]

    objectNameString = (" ".join(objectNameReadLine.split())).split(delimiter)
    objectNameString = objectNameString[1]

    RAString = (" ".join(RAReadLine.split())).split(delimiter)
    RAString = RAString[1]

    DECString = (" ".join(DECReadLine.split())).split(delimiter)
    DECString = DECString[1]

    dataObservedByString = (" ".join(dataObservedByReadLine.split())).split(delimiter)
    dataObservedByString = dataObservedByString[1]

    epochTimeMJDString = (" ".join(epochTimeMJDReadLine.split())).split(delimiter)
    epochTimeMJDString = epochTimeMJDString[1]

    barycenterFlagString = (" ".join(barycenterFlagReadLine.split())).split(delimiter)
    barycenterFlagString = barycenterFlagString[1]

    nbinsTimeSeriesString = (" ".join(nbinsTimeSeriesReadLine.split())).split(delimiter)
    nbinsTimeSeriesString = nbinsTimeSeriesString[1]

    timeResolutionString = (" ".join(timeResolutionReadLine.split())).split(delimiter)
    timeResolutionString = timeResolutionString[1]

    dataBreakFlagString = (" ".join(dataBreakFlagReadLine.split())).split(delimiter)
    dataBreakFlagString = dataBreakFlagString[1]

    EMBandString = (" ".join(EMBandReadLine.split())).split(delimiter)
    EMBandString = EMBandString[1]

    beamDiameterString = (" ".join(beamDiameterReadLine.split())).split(delimiter)
    beamDiameterString = beamDiameterString[1]

    dispersionMeasureString = (" ".join(dispersionMeasureReadLine.split())).split(delimiter)
    dispersionMeasureString = dispersionMeasureString[1]

    centralFreqBottomBandString = (" ".join(centralFreqBottomBandReadLine.split())).split(delimiter)
    centralFreqBottomBandString = centralFreqBottomBandString[1]

    bandwidthString = (" ".join(bandwidthReadLine.split())).split(delimiter)
    bandwidthString = bandwidthString[1]

    numChannelsString = (" ".join(numChannelsReadLine.split())).split(delimiter)
    numChannelsString = numChannelsString[1]

    channelBandwidthString = (" ".join(channelBandwidthReadLine.split())).split(delimiter)
    channelBandwidthString = channelBandwidthString[1]

    dataAnalyzedByString = (" ".join(dataAnalyzedByReadLine.split())).split(delimiter)
    dataAnalyzedByString = dataAnalyzedByString[1]

    readFile.close()

#    print(dataFileNameString)
#    print(telescopeNameString)
#    print(instrumentString)
#    print(objectNameString)
#    print(RAString)
#    print(DECString)
#    print(dataObservedByString)
#    print(epochTimeMJDString)
#    print(barycenterFlagString)
#    print(nbinsTimeSeriesString)
#    print(timeResolutionString)
#    print(dataBreakFlagString)
#    print(EMBandString)
#    print(beamDiameterString)
#    print(dispersionMeasureString)
#    print(centralFreqBottomBandString)
#    print(bandwidthString)
#    print(numChannelsString)
#    print(channelBandwidthString)
#    print(dataAnalyzedByString)


    infFileData = infDataStruct(dataFileNameString, telescopeNameString, \
    instrumentString, objectNameString, RAString, DECString, dataObservedByString, \
    epochTimeMJDString, barycenterFlagString, nbinsTimeSeriesString, \
    timeResolutionString, dataBreakFlagString, EMBandString, beamDiameterString, \
    dispersionMeasureString, centralFreqBottomBandString, bandwidthString, \
    numChannelsString, channelBandwidthString, dataAnalyzedByString)

    return infFileData;



""" Write the data from PRESTO's .dat file and the information from the .inf
    file to a new filterbank file. """
def convertPrestoDatFileToFilterbank(prestoData, infFileData, filterbankFilename,
                                     prestoDatFilename):

    fb = fil.Filterbank(filterbankFilename, prestoData)


    rawdatafile = os.path.split(filterbankFilename)[1]
    source_name = infFileData.objectNameString
    telescope_id = infFileData.telescopeNameString

    if (telescope_id == "Canberra"):
        telescope_id = 12
    else:
        telescope_id = 63

    machine_id = 999
    data_type = 1
    nifs = 1
    tsamp = float(infFileData.timeResolutionString)
    barycentric = int(infFileData.barycenterFlagString)

    fb.rawdatafile = rawdatafile
    fb.source_name = source_name
    fb.telescope_id = telescope_id
    fb.machine_id = machine_id
    fb.data_type = data_type
    fb.nifs = nifs
    fb.tsamp = tsamp
    fb.barycentric = barycentric
    fb.setRaDec(infFileData.RAString, infFileData.DECString)
    fb.tstart = float(infFileData.epochTimeMJDString)
    fb.nbeams = 0.0
    fb.nchans = 1
    fb.foff =  -float(infFileData.bandwidthString)

    fch1 = np.divide(np.add(np.add(float(infFileData.centralFreqBottomBandString),
                                   np.subtract(float(infFileData.numChannelsString), 1.0)),
                            float(infFileData.centralFreqBottomBandString)),
                     2.0)
    fb.fch1 = fch1
    fb.nbits = 32

    fb.writeHeader(filterbankFilename)

#    call("filedit --ra %s %s" % (infFileData.RAString.replace(":", ""), filterbankFilename), shell=True)
#    call("filedit --dec %s %s" % (infFileData.DECString.replace(":", ""), filterbankFilename), shell=True)
#    call("filedit --tstart %s %s" % (infFileData.epochTimeMJDString, filterbankFilename), shell=True)
#    call("filedit --beam 0 %s" % filterbankFilename, shell=True)
#    call("filedit --nbeams 0 %s" % filterbankFilename, shell=True)
#    call("filedit --nchan 1 %s" % filterbankFilename, shell=True)
#
#    fch1 = np.divide(np.add(np.add(float(infFileData.centralFreqBottomBandString),
#                            np.subtract(float(infFileData.numChannelsString), 1.0)),
#                            float(infFileData.centralFreqBottomBandString)), 2.0)
#
#    call("filedit --fch1 %.8f %s" % (fch1, filterbankFilename), shell=True)
#
#    call("filedit --foff %.10f %s" % (fb.foff, filterbankFilename), shell=True)
#    call("filedit --nbits 32 %s" % filterbankFilename, shell=True)


    fb.rawData = prestoData

    fb.writeData(filterbankFilename)

#    x = fb.readDataChunk(22608530)
#    x = fb.readDataChunk(10000)
#    plt.plot(x)
#    plt.savefig("test.png")



    return;












def main():
    if (len(sys.argv) != 4):
        print("##################################")
        print("Aaron B. Pearlman")
        print("aaron.b.pearlman@caltech.edu")
        print("Division of Physics, Mathematics, and Astronomy")
        print("California Institute of Technology")
        print("Jet Propulsion Laboratory")
        print("##################################\n")

        print("m_convertPRESTODatFileToFilterbank.py - Read in a PRESTO .dat " \
        "file and the information from its corresponding .inf file. Convert the " \
        "data into filterbank format.\n")

        print("Usage: m_convertPRESTODatFileToFilterbank.py [Input PRESTO .dat filename] [Input PRESTO .inf filename] [Output filterbank filename]\n")
        print("Example: m_convertPRESTODatFileToFilterbank.py slcp_16a245_b1.0_rfi1s.dat slcp_16a245_b1.0_rfi1s.inf slcp_16a245_b1.0_rfi1s_PRESTO.corr")
        sys.exit(0)

    prestoDatFilename = sys.argv[1]
    prestoInfFilename = sys.argv[2]
    filterbankFilename = sys.argv[3]



    # Read in data from the PRESTO .dat file.
    prestoSamples, prestoData = readPrestoDatFile(prestoDatFilename)

    # Read in information from the PRESTO .inf file.
    infFileData = readPrestoInfFile(prestoInfFilename)

    # Convert the PRESTO .dat file to a filterbank file.
    convertPrestoDatFileToFilterbank(prestoData, infFileData, filterbankFilename,
                                     prestoDatFilename)



if __name__ == "__main__":
    main()
