#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 16:52:32 2017

@author: Aaron B. Pearlman
         aaron.b.pearlman@caltech.edu
         Division of Physics, Mathematics, and Astronomy
         California Institute of Technology
         Jet Propulsion Laboratory

m_filter_fft_rect_window.py - Input is a non-barycentered filterbank file.

                              The filterbank file is split into single channel
                              .dat files using fb_truncate.py and PRESTO's
                              prepdata tool. Each .dat file is FFT'd using
                              PRESTO's realfft command and an input zap mask is
                              applied to filter out undesired signals and
                              harmonics in the Fourier domain using a
                              rectangular window function. After these
                              undesired signals are "zapped", the FFT'd data is
                              then inverse FFT'd and the resultant time-series
                              (now overwritten in the .dat file) are converted
                              into single channel filterbank files.

                              A filtered version of the input non-barycentered
                              filterbank file using these processing steps is
                              output at the end of processing using a custom
                              version of SIGPROC's splice command.

                              NOTE: The data are not barycentered and we use
                              PRESTO's prepdata command with -nobary, -noclip,
                              and DM=0.0 are used.

"""

import numpy as np
import sys
import os
from subprocess import call, check_output

import filterbank_class as fil




""" Create single channel barycentered filterbank files using PRESTO. """
def filterbankBarySingleChannels(filterbankFilename, zapListFilename,
                                 topFrequency, bottomFrequency, offsetFrequency,
                                 rootOutputFilename, DM=0.0):

    print("\n\nCreating barycentered single channel filterbank files...\n\n")

    progressCounter = 1.0
    progressTotal = len(np.arange(topFrequency, bottomFrequency + offsetFrequency, offsetFrequency))

    barycenteredFilterbankFilenames = []

    for frequency in np.arange(topFrequency, bottomFrequency + offsetFrequency, offsetFrequency):


        print("\n\n[Barycentering Single Channel Filterbank Files]: %.2f%% Completed\n\n"
              % np.multiply(np.divide(progressCounter, progressTotal), 100.0))


        singleChannelOutputFilterbankFilename = "%s_%iMHz.corr" % (rootOutputFilename, 10.*frequency)

        # Run fb_truncate.py to split the filterbank file into single channels.
        call("fb_truncate.py %s -H %f -L %f -o %s" \
             % (filterbankFilename, frequency, frequency, singleChannelOutputFilterbankFilename), shell=True)



        prepdataOutputFilename = os.path.splitext(singleChannelOutputFilterbankFilename)[0]

        # Run PRESTO to generate barycentered time series for each single channel
        # filterbank file.
        if (DM == 0.0):
            call("prepdata -nobary -noclip -filterbank %s -o %s" % (singleChannelOutputFilterbankFilename, prepdataOutputFilename), shell=True)
            call("realfft %s.dat" % prepdataOutputFilename, shell=True)
            call("zapbirds -zap -zapfile %s %s.fft" % (zapListFilename, prepdataOutputFilename), shell=True)
            call("realfft -inv %s.fft" % prepdataOutputFilename, shell=True)

        else:

            # In practice, this is not used because we build the barycenter
            # filterbank file using DM = 0.0 and then dedisperse later.
            call("prepdata -noclip -dm %.2f -filterbank %s -o %s" % (DM, singleChannelOutputFilterbankFilename, prepdataOutputFilename),
                 shell=True)



        # Clean up the single channel filterbank files.
        call("rm %s" % singleChannelOutputFilterbankFilename, shell=True)



        prepdataOutputDatFilename = prepdataOutputFilename + ".dat"
        prepdataOutputInfFilename = prepdataOutputFilename + ".inf"
        prepdataOutputFilterbankBaryFilename = prepdataOutputFilename + "_bary.corr"

        # Convert the barycentered PRESTO .dat files into barycentered
        # single channel filterbank files.
        call("m_convertPRESTODatFileToFilterbank.py %s %s %s" \
             % (prepdataOutputDatFilename, prepdataOutputInfFilename,
                prepdataOutputFilterbankBaryFilename), shell=True)

        # Clean up the PRESTO .dat and .inf files and .fft.
        call("rm %s" % prepdataOutputDatFilename, shell=True)
        call("rm %s" % prepdataOutputInfFilename, shell=True)
        call("rm %s" % (prepdataOutputFilename + ".fft"), shell=True)


        barycenteredFilterbankFilenames.append(prepdataOutputFilterbankBaryFilename)

        # Increamenet the progress counter.
        progressCounter = progressCounter + 1

    return barycenteredFilterbankFilenames;



""" Merge the barycentered single channel filterbank files from PRESTO back
    together. """
def mergeBaryFilterbankFiles(barycenteredFilterbankFilenames, rootOutputFilename):

    spliceBlockSize = 500
    spliceBlocks = np.linspace(0, len(barycenteredFilterbankFilenames), len(barycenteredFilterbankFilenames) // spliceBlockSize + 2, dtype=int)

    spliceOutputBaryFilterbankFilename = []
    spliceFileListString = ""
    for spliceIndex in np.arange(0, len(spliceBlocks) - 1, 1):

        fileListString = ""

        print("\n\nSplicing the barycentered single channel filterbank files together into a single filterbank file...\n\n")

        for fileListIndex in np.arange(spliceBlocks[spliceIndex], spliceBlocks[spliceIndex + 1], 1):

            fileListString = fileListString + barycenteredFilterbankFilenames[fileListIndex] + " "

        outputBaryFilterbankFilename = rootOutputFilename + "_" + str(spliceIndex) + "_filter.corr"

        # Combine all of the barycentered single channel filterbank files from
        # PRESTO into a single barycentered filterbank file using sigproc's splice
        # tool (customized).
#        call("splice %s > %s" % (fileListString, outputBaryFilterbankFilename),
#             shell=True)
        call("splice_fb_bary %s-o %s" % (fileListString, outputBaryFilterbankFilename),
             shell=True)


        # Clean up the barycentered single channel filterbank files from PRESTO.
        for fileListIndex in np.arange(spliceBlocks[spliceIndex], spliceBlocks[spliceIndex + 1], 1):
            call("rm %s" % barycenteredFilterbankFilenames[fileListIndex], shell=True)

        spliceOutputBaryFilterbankFilename.append(outputBaryFilterbankFilename)
        spliceFileListString = spliceFileListString + outputBaryFilterbankFilename + " "


    outputBaryFilterbankFilename = rootOutputFilename + "_filter.corr"

    # Combine all of the barycentered single channel filterbank files from
    # PRESTO into a single barycentered filterbank file using sigproc's splice
    # tool (customized).
#    call("splice %s > %s" % (spliceFileListString, outputBaryFilterbankFilename),
#         shell=True)
    call("splice_fb_bary %s-o %s" % (spliceFileListString, outputBaryFilterbankFilename),
         shell=True)

    for fileListIndex in np.arange(0, len(spliceOutputBaryFilterbankFilename), 1):
        call("rm %s" % spliceOutputBaryFilterbankFilename[fileListIndex], shell=True)

    return;



def main():
    if (len(sys.argv) != 3):
        print("##################################")
        print("Aaron B. Pearlman")
        print("aaron.b.pearlman@caltech.edu")
        print("Division of Physics, Mathematics, and Astronomy")
        print("California Institute of Technology")
        print("Jet Propulsion Laboratory")
        print("##################################\n")

        print("m_filter_fft_rect_window.py - Input is a non-barycentered " \
        "filterbank file.\n\n" \
        " The filterbank file is split into single channel .dat files using " \
        "fb_truncate.py and PRESTO's prepdata tool. Each .dat file is FFT'd " \
        "using PRESTO's realfft command and an input zap mask is applied to " \
        "filter out undesired signals and harmonics in the Fourier domain " \
        "using a rectangular window function. After these undesired signals " \
        "are \"zapped\", the FFT'd data is then inverse FFT'd and the resultant " \
        "time-series (now overwritten in the .dat file) are converted into " \
        "single channel filterbank files.\n\n" \

        "A filtered version of the input non-barycentered filterbank file " \
        "using these processing steps is output at the end of processing " \
        "using a custom version of SIGPROC's splice command.\n\n" \

        "NOTE: The data are not barycentered and we use PRESTO's prepdata " \
        "command with -nobary, -noclip, and DM=0.0.\n\n" \

        "NOTE2 (**): You first need to run PRESTO's makezaplist.py sband.birds, " \
        "where sband.birds is a sample *.birds file. A sband.inf file also " \
        "needs to be created.\n")

        print("Usage: m_filter_fft_rect_window.py [Input filterbank filename] [Input zapfile filename]\n")
        print("Example: m_filter_fft_rect_window.py srcp.corr srcp.zaplist")
        sys.exit(0)

    filterbankFilename = sys.argv[1]
    zapListFilename = sys.argv[2]


    # Split the filterbank into barycentered single channel filterbank files.
    print("Reading header from raw filterbank file...\n\n")
    fb = fil.Filterbank(filterbankFilename, readAllFlag=0)

    topFrequency = float(fb.fch1)
    offsetFrequency = fb.foff
    bandwidth = float(fb.nchans)*offsetFrequency
    bottomFrequency = topFrequency + (bandwidth - offsetFrequency)
    rootOutputFilename = os.path.splitext(filterbankFilename)[0]


    # Default: Barycenter using PRESTO's prepdata command with DM = 0.0.
    barycenteredFilterbankFilenames = filterbankBarySingleChannels(filterbankFilename, zapListFilename,
                                                                   topFrequency, bottomFrequency, offsetFrequency,
                                                                   rootOutputFilename)

    mergeBaryFilterbankFiles(barycenteredFilterbankFilenames, rootOutputFilename)



if __name__ == "__main__":
    main()
