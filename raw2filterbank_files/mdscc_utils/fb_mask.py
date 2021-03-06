#!/usr/bin/env python

"""
A script to mask a filterbank file in frequency.

Shahab Arabshahi 
based on fb_truncate.py script by Patrick Lazarus, June 5, 2017
"""

import sys
import copy
import os.path
import warnings
import optparse

import numpy as np

import filterbank

BLOCKSIZE = 1e5 # Number of spectra to manipulate at once

def main():
    infn = args[0]
    print "Reading filterbank file (%s)" % infn
    fil = filterbank.FilterbankFile(infn)
    startbin = 0
    endbin = fil.nspec
    nspec = endbin-startbin
   
    # Determine lo/hi channels to mask
    # If high frequencies come first in spectra 'hichan' refers to 
    # the lo-freq cutoff and 'lochan' refers to the hi-freq cutoff.
    if options.lo_freq is None:        
        raise ValueError("Low frequency is needed.")
    else:
        ichan = int(np.round((options.lo_freq-fil.fch1)/fil.foff))
#        if ichan < 0 or ichan >= fil.nchans:
#            raise ValueError("Invalid Low frequency!")
        if fil.foff > 0:
            lochan = ichan
        else:
            hichan = ichan+1
    if options.hi_freq is None:
        raise ValueError("High frequency is needed.")
    else:
        ichan = int(np.round((options.hi_freq-fil.fch1)/fil.foff))
#        if ichan < 0 or ichan >= fil.nchans:
#            raise ValueError("Invalid High frequency!")
        if fil.foff > 0:
            hichan = ichan+1
        else:
            lochan = ichan
            
    if lochan < 0:
        lochan = 0
    if hichan >= fil.nchans:
        hichan = fil.nchans - 1


    print "Will extract"
    print "    %d bins (%d to %d incl.)" % (nspec, startbin, endbin-1)
    print "    (Original num bins: %d)" % fil.nspec
    print "    (Original num chans: %d)" % fil.nchans

    # Create output file
    outfn = options.outname % fil.header
    print "Creating out file: %s" % outfn
    outhdr = copy.deepcopy(fil.header)
    filterbank.create_filterbank_file(outfn, outhdr, nbits=fil.nbits)
    outfil = filterbank.FilterbankFile(outfn, mode='write')

    # Write data
    sys.stdout.write(" %3.0f %%\r" % 0)
    sys.stdout.flush()
    nblocks = int(nspec/options.block_size)
    remainder = nspec % options.block_size
    oldprogress = -1
    for iblock in np.arange(nblocks):
        lobin = iblock*options.block_size + startbin
        hibin = lobin+options.block_size
        spectra = fil.get_spectra(lobin, hibin)
        spectra[:,lochan:hichan] = 0.0 # zero channels
        outfil.append_spectra(spectra)
        progress = int(100.0*((hibin-startbin)/nspec))
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    # Read all remaining spectra
    if remainder:
        spectra = fil.get_spectra(endbin-remainder, endbin)
        spectra[:,lochan:hichan] = 0.0 # zero channels
        outfil.append_spectra(spectra)
    sys.stdout.write("Done   \n")
    sys.stdout.flush()


if __name__ == '__main__':
    parser = optparse.OptionParser(prog='fb_truncate.py', \
                    version="v0.1 Patrick Lazarus (Aug. 28, 2012)")
    parser.add_option("-L", "--lo-freq", dest="lo_freq", type='float', \
                    help="Desired low frequency for output file. Note: " \
                        "actual low frequency will be rounded to the nearest" \
                        "channel (Default: Don't truncate low-freq channels)", \
                    default=None)
    parser.add_option("-H", "--hi-freq", dest="hi_freq", type='float', \
                    help="Desired high frequency for output file. Note: " \
                        "actual high frequency will be rounded to the nearest" \
                        "channel (Default: Don't truncate high-freq channels)", \
                    default=None)
    parser.add_option("--block-size", dest='block_size', default=BLOCKSIZE, \
                    type='float', \
                    help="Number of spectra per block. This is the amount " \
                        "of data manipulated/written at a time. (Default: " \
                        " %d spectra)" % BLOCKSIZE)
    parser.add_option("-o", "--outname", dest='outname', action='store', \
                    help="The name of the output file.")
    (options, args) = parser.parse_args()
    if not hasattr(options, 'outname'):
        raise ValueError("An output file name _must_ be provided!")
    main()
