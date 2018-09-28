import os, time
from mdscc_raw2filterbank_class import *

##### Calculate SK Estimator scale factor #####
def process_SK_ScaleCalc(filename, header, M):

    print "Start calculating SK estimator scale factor"

    NCH = header['nchans']

    S1ibuf = FourChannelData()
    S2ibuf = FourChannelData()


    print "filename = ", filename
    statinfo = os.stat(filename)
    if statinfo.st_size == 0:
        raise IOError("Zero file size!")

    f_in= open(filename)


    ssize = (NCH * 4 * 2) * 2   # ( Number of freq channels x
                                # number of chs x number of moments )
                                # x number of bytes per data

    fsize = statinfo.st_size
    nspectra = fsize / ssize - 32

    S1ibuf.initialize(NCH, nspectra)
    S2ibuf.initialize(NCH, nspectra)


    stime = time.time()
    print "%d spectra to read ..." %nspectra
    ##### Read raw data block:
    read_4chan_spectra(S1ibuf,S2ibuf,f_in, nspectra)
    print "read file in %f seconds ..." %(time.time()-stime)

    ##### Flip 1d data block:
    if header['foff'] < 0:
        S1ibuf.flip()
        S2ibuf.flip()


    f_in.close()


    f_scale =  FourChannelData_SK_CalcScale(S1ibuf,S2ibuf,M)

    print "End calculating SK estimator scale factor"

    return f_scale
###############################################


def process_raw2filterbank(infiles, outpath, header, f_scale, FCH1, M, output_flags):


    NCH = header['nchans']

    headers = FourChannelHeader(header, FCH1)

    S1 = FourChannelData()
    S2 = FourChannelData()
    SK = FourChannelData()

    ##### Output raw data as filterbank files #####
    if output_flags['S1raw'] :
        S1.initOutput(outpath, header['bname'], 'S1')
        # write header
        S1.writeHeader(headers)

    if output_flags['S2raw'] :
        S2.initOutput(outpath, header['bname'], 'S2')
        # write header
        S2.writeHeader(headers)
    ###############################################

    ##### Output SK Estimator values as filterbank files #####
    if output_flags['SKestimator']:
        SK.initOutput(outpath, header['bname'], 'SK')
        # write header
        SK.writeHeader(headers)
    ##########################################################

    for filename in infiles:

        print "filename = ", filename
        statinfo = os.stat(filename)
        if statinfo.st_size == 0:
            raise IOError("Zero file size!")


        f_in= open(filename)


        ssize = (NCH * 4 * 2) * 2   # ( Number of freq channels x
                                    # number of chs x number of moments )
                                    # x number of bytes per data

        fsize = statinfo.st_size
        nspectra = fsize / ssize - 32

        S1.initialize(NCH, nspectra)
        S2.initialize(NCH, nspectra)
        SK.initialize(NCH, nspectra)


        stime = time.time()
        print "%d spectra to read ..." %nspectra
        read_4chan_spectra(S1,S2,f_in,nspectra)
        print "read file in %f seconds ..." %(time.time()-stime)

        ##### Flip 1d data block:
        if header['foff'] < 0:
            S1.flip()
            S2.flip()


        ##### Output raw data as filterbank files #####
        if output_flags['S1raw']:
            stime = time.time()
            S1.writeOutput()
            print "wrote S1 raw data block wrote in %s seconds" %(time.time()-stime)

        if output_flags['S2raw']:
            stime = time.time()
            S2.writeOutput()
            print "wrote S2 raw data block wrote in %s seconds" %(time.time()-stime)
        ###############################################


        ##### Output SK Estimator values as filterbank files #####
        if output_flags['SKestimator']:
            ##### Calculate SK Estimator:
            stime = time.time()
            SK.SKCalc(S1, S2, f_scale, M)
            print "calculated SK Estimator in %s seconds" %(time.time()-stime)

            stime = time.time()
            SK.writeOutput()
            print "wrote SK Estimator data block in %s seconds" %(time.time()-stime)
        ##########################################################


        f_in.close()

    ##### Close raw data filterbank files #####
    if output_flags['S1raw']:
        S1.finishOutput()

    if output_flags['S2raw']:
        S2.finishOutput()
    ###############################################

    ##### Close SK Estimator values filterbank files #####
    if output_flags['SKestimator']:
        SK.finishOutput()
    ##########################################################