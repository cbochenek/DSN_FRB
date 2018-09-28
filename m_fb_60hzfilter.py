#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 19:14:37 2016

@author: pearlman

Aaron B. Pearlman
aaron.b.pearlman@caltech.edu
Division of Physics, Mathematics, and Astronomy
California Institute of Technology
Jet Propoulsion Laboratory

m_fb_60hzfilter.py - Script for filtering out 60 hz instrumental noise and
                     higher frequency harmonics.

"""

import getopt, sys
from subprocess import call, check_call, check_output, Popen

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import h5py

sys.path.append("/Users/aaron/pulsar_software/presto/lib/python/")
import filterbank

from scipy import signal



BLOCKSIZE = 1e6

""" Read the filterbank file into memory. Store the data in a dynamically
    accessible h5py file, stored in a binary .hdf5 file. """
def readFilterbank(inputFilename, logFile=""):

    if (logFile == ""):
        print("Reading filterbank file (%s)...\n" % inputFilename)
    else:
        logFile.write("Reading filterbank file (%s)...\n\n" % inputFilename)

    fb = filterbank.FilterbankFile(inputFilename)

    inputHeader = copy.deepcopy(fb.header)
    inputNbits = fb.nbits

    totalChans = fb.nchans
    nchans = np.arange(0, fb.nchans-1, 1) # Top of the band is index 0.

    freqs = fb.frequencies

    startbin = 0
    endbin = fb.nspec
    nspec = np.subtract(endbin, startbin)

    nblocks = int(np.divide(nspec, BLOCKSIZE))
    remainder = nspec % BLOCKSIZE
    totalBlocks = nblocks

    if (remainder):
        totalBlocks = nblocks + 1



    h5pyFile = h5py.File("%s.hdf5" % inputFilename, "w")
    spectraData = h5pyFile.create_dataset("data", (totalChans, nspec), dtype="float32")

    
    iblock = 0
    
    for iblock in np.arange(0, nblocks, 1):

        progress = np.multiply(np.divide(iblock + 1.0, totalBlocks), 100.0)

        if (logFile == ""):
            sys.stdout.write("Reading... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Reading... [%3.2f%%]\n" % progress)


        lobin = int(np.add(np.multiply(iblock, BLOCKSIZE), startbin))
        hibin = int(np.add(lobin, BLOCKSIZE))

        spectra = fb.get_spectra(lobin, hibin)


        for ichan in np.arange(0, totalChans, 1):

            spectraData[ichan, lobin:hibin] = spectra[:,ichan]


    if (remainder):

        progress = np.multiply(np.divide(iblock + 2.0, totalBlocks), 100.0)

        if (logFile == ""):
            sys.stdout.write("Reading... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Reading... [%3.2f%%]\n" % progress)


        lobin = int(np.subtract(endbin, remainder))
        hibin = int(endbin)

        spectra = fb.get_spectra(lobin, hibin)

        for ichan in np.arange(0, totalChans, 1):

            spectraData[ichan, lobin:hibin] = spectra[:,ichan]

    if (logFile == ""):
        print("\n")
    else:
        logFile.write("\n")

    return spectraData, inputHeader, inputNbits, h5pyFile;



""" Write the filterbank data from memory to a filterbank file. """
def writeFilterbank(outputFilename, spectraData, inputHeader, inputNbits,
                    logFile=""):

    if (logFile == ""):
        print("Writing filterbank file (%s)...\n" % outputFilename)
    else:
        logFile.write("Writing filterbank file (%s)...\n\n" % outputFilename)


    filterbank.create_filterbank_file(outputFilename, inputHeader, nbits=inputNbits)
    outfil = filterbank.FilterbankFile(outputFilename, mode='write')

    startbin = 0
    endbin = np.shape(spectraData)[1]

    nblocks = int(np.divide(endbin, BLOCKSIZE))
    remainder = endbin % BLOCKSIZE
    totalBlocks = nblocks

    if (remainder):
        totalBlocks = nblocks + 1


    
    iblock = 0
    
    for iblock in np.arange(0, nblocks, 1):

        progress = np.multiply(np.divide(iblock + 1.0, totalBlocks), 100.0)

        if (logFile == ""):
            sys.stdout.write("Writing... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Writing... [%3.2f%%]\n" % progress)



        lobin = int(np.add(np.multiply(iblock, BLOCKSIZE), startbin))
        hibin = int(np.add(lobin, BLOCKSIZE))

        spectra = spectraData[:,lobin:hibin].T
        outfil.append_spectra(spectra)


    if (remainder):

        progress = np.multiply(np.divide(iblock + 2.0, totalBlocks), 100.0)

        if (logFile == ""):
            sys.stdout.write("Writing... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Writing... [%3.2f%%]\n" % progress)



        lobin = int(np.subtract(endbin, remainder))
        hibin = int(endbin)

        spectra = spectraData[:,lobin:hibin].T
        outfil.append_spectra(spectra)

    if (logFile == ""):
        print("\n")
    else:
        logFile.write("\n")

    return;


""" Create a butterworth bandstop filter (use an analog filter for real data!). """
def butter_bandstop(nyq, cutoff_freq_start, cutoff_freq_stop, order=3):
    
    cutoff_freq_start = cutoff_freq_start / nyq
    cutoff_freq_stop = cutoff_freq_stop / nyq
    
    b, a = signal.butter(order, [cutoff_freq_start, cutoff_freq_stop], btype="bandstop", analog=False)
    
    return b, a



""" Filter 60 Hz signal and higher order harmonics from the fb file.
    Typical attenuation is ~200 dB around the filtered frequencies. """
def fb_filter_60Hz(fb_data, fb_header, maxHarmonicFrequency):
    
    timeRes = float(fb_header["tsamp"])
    nsamples = np.shape(fb_data)[1]
    duration = np.divide(np.multiply(nsamples, timeRes), 3600.0)
    
    print("Time Resolution: %.6f s" % timeRes)
    print("nsamples: %i" % nsamples)
    print("Duration: %.2f hr\n" % duration)
    
    
    fs = 1.0 / timeRes
    nyq = 0.5 * fs
    
    cutoff_freq_fundamental = 60.0
    cutoff_freq_start = 0.0
    cutoff_freq_stop = 0.0
    
    for filterFreq in np.arange(cutoff_freq_fundamental, np.add(maxHarmonicFrequency, cutoff_freq_fundamental), cutoff_freq_fundamental):
        
        print("Filtering frequency: %.1f Hz" % filterFreq)
        
        if (filterFreq == cutoff_freq_fundamental):
            
            cutoff_freq_start = 58.
            cutoff_freq_stop = 62.
        
        else:
            
            cutoff_freq_start = np.subtract(filterFreq, 2.)
            cutoff_freq_stop = np.add(filterFreq, 2.)
        
        b, a = butter_bandstop(nyq, cutoff_freq_start, cutoff_freq_stop, order=5)
        w, h = signal.freqz(b, a, worN=100000)
        
	for index in np.arange(0, np.shape(fb_data)[0], 1):
        
            fb_data[index] = signal.filtfilt(b, a, fb_data[index])

    return fb_data;



def usage():
    print("##################################")
    print("Aaron B. Pearlman")
    print("aaron.b.pearlman@caltech.edu")
    print("Division of Physics, Mathematics, and Astronomy")
    print("California Institute of Technology")
    print("Jet Propulsion Laboratory")
    print("##################################\n")

    print """
    usage:  m_fb_60hzfilter.py [options]
        [-h, --help]                    : Display this help
        [--inputFilename]               : Name of input filterbank file
        [--outputFilename]              : Name of output filterbank file created after
                                          filtering is completed
        [--maxHarmonicFrequency]        : Filter harmonics up to maxHarmonicFrequency (Hz)
        [--outputDir]                   : Output directory where the products of the data
                                          analysis will be stored
        [--clean]                       : Flag to clean up intermediate reduction products.
                                          Default is FALSE
        
        
        This program reads a filterbank file, stores it in dynamic memory,
        and then removes 60 Hz instrumental noise from each channel of the
        filterbank file. Additional harmonics of this noise signal, up to
        [maxHarmonicFrequency] Hz are also filtered. The filtered data is
        written to a new filterbank file.
        
        TO DO: A parallelized version of this algorithm is currently under
        development (will be completed circa Sept. 2018).
        
    Example: m_fb_60hzfilter.py --inputFilename input.corr --outputFilename output.corr --maxHarmonicFrequency 200.0 --outputDir /home/pearlman/fb_data/ --clean True

    """





def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "inputFilename:outputFilename:timeConstLong:timeConstShort:numProcessors:outputDir:logFile:clean:",
                                   ["help", "inputFilename=", "outputFilename=",
                                    "maxHarmonicFrequency=", "outputDir=",
                                    "clean"])
    
    except getopt.GetoptError:
        # Print help information and exit.
        usage()
        sys.exit(2)
    
    if (len(sys.argv) == 1):
        usage()
        sys.exit(2)
    
    inputFilename=None
    outputFilename=None
    maxHarmonicFrequency=None
    outputDir=None
    clean=None
    
    for o, a in opts:
        if (o in ("-h", "--help")):
            usage()
            sys.exit()
        
        if o in ("--inputFilename"):
            inputFilename = a
        if o in ("--outputFilename"):
            outputFilename = a
        if o in ("--maxHarmonicFrequency"):
            maxHarmonicFrequency = a
        if o in ("--outputDir"):
            outputDir = a
        if o in ("--clean"):
            clean = True
    
    if ((inputFilename == None) | (outputFilename == None) \
        | (maxHarmonicFrequency == None)):
        usage()
        sys.exit()
    
    if (maxHarmonicFrequency != None):
        maxHarmonicFrequency = float(maxHarmonicFrequency)
    
    
    
    # Read the filterbank data.
    fb_data, fb_header, inputNbits, h5pyFile = readFilterbank(inputFilename)
    
    # Filter out the 60 Hz noise and higher frequency harmonics.
    fb_data = fb_filter_60Hz(fb_data, fb_header, maxHarmonicFrequency)
    
    # Write the filtered filterbank data to a new file.
    writeFilterbank(outputFilename, fb_data, fb_header, inputNbits,
                    logFile="")
    
    
    
    
    if (clean == True):
        
        if (outputDir != None):
            call("rm -rf %s/*hdf5" % outputDir, shell=True)
        else:
            call("rm -rf *hdf5", shell=True)
    
    
    
if __name__ == "__main__":
    main()
