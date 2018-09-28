#! /Users/majid/anaconda/bin/python

#
# example:
#  fn = "/Users/majid/psr_recording_000_201409121751"
#  import dsn_read_incoh
#  ps,ks = dsn_read_incoh.read_incoh(fn)
#  dsn_read_incoh.ds_diag(ps)
#

# read incoherent mode data

import numpy as np
from itertools import chain
from collections import deque
from itertools import islice
import sys
import struct
import os
import glob
import time
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy.fftpack import fft
import os.path
import getopt

import proc_params


# set default header info here
#telescope_id = 999
telescope_id = 12
machine_id = 999

data_type = 1 # filterbank
#foff = -1.0
#nchans = 1024
obits = 32

nifs = 1
barycentric = 0

def main(idir,odir, source, tsamp, tag, nproc):
    # source   = proc_params.SOURCE
    tsamp    = tsamp * 1.0e-6 # convert usec to sec
    #tsamp    = (proc_params.TSAMP) * 1.0e-6 # convert usec to sec
    ch1freq  = proc_params.CH1FREQ
    ch2freq  = proc_params.CH2FREQ
    ch3freq  = proc_params.CH3FREQ
    ch4freq  = proc_params.CH4FREQ
    ch1label = proc_params.CH1LABEL
    ch2label = proc_params.CH2LABEL
    ch3label = proc_params.CH3LABEL
    ch4label = proc_params.CH4LABEL
    bw       = proc_params.BW
    nskip    = proc_params.NSKIP

    filePattern = idir + '/*.dada'
    files = sorted(glob.iglob(filePattern))
    first = 1

    nfiles_processed = 0

    for fn in files:
        print 'processing %s ...' % (fn)
        if first:
            # write headers
            mjd, nsamps, bname = parse_dada_name(fn)
            print("main: mjd %f nsamps %d bname %s" %(mjd, nsamps, bname))

            bname = tag + "-" + bname

            # set the names of the filterbank files using the first dada filename
            ofn1 = odir + "/" + bname + "_" + ch1label + ".fil"
            ofn2 = odir + "/" + bname + "_" + ch2label + ".fil"
            ofn3 = odir + "/" + bname + "_" + ch3label + ".fil"
            ofn4 = odir + "/" + bname + "_" + ch4label + ".fil"

            # BW = 1000 # MHz ==> channel bw is then 976562.5 Hz or slightly less than 1 MHz
            write_header(ofn1, ch1freq+bw, source, tsamp, mjd)
            write_header(ofn2, ch2freq+bw, source, tsamp, mjd)
            write_header(ofn3, ch3freq+bw, source, tsamp, mjd)
            write_header(ofn4, ch4freq+bw, source, tsamp, mjd)

        # split into 4 files
        dada_to_4ch(fn, ofn1, ofn2, ofn3, ofn4)
        nfiles_processed = nfiles_processed + 1
        first = 0
        if (nproc != 0) and (nproc == nfiles_processed):  # process only NPROC files
            break


def dada_to_4ch(ifn, ch1,ch2,ch3,ch4):
    t0 = time.time()
    UTC_OFFSET_TIMEDELTA = dt.datetime.utcnow() - dt.datetime.now()
    firstSpectra = 1
    myname = "dada_4ch"
    if (proc_params.DEBUG):
        print("%s: ifn: %s" %(myname, ifn))

    ofp1 = open(ch1, "ab")
    ofp2 = open(ch2, "ab")
    ofp3 = open(ch3, "ab")
    ofp4 = open(ch4, "ab")

    nspec = get_nspec_dada(ifn)
    print("%s: nspec: %d" %(myname, nspec))

    ifp, s1, s2, s3, s4 = open_file(ifn)

    foff     = proc_params.FOFF

    for i in range(nspec):

        #if (i % 20000) == 0:
        #    print "%s: spectra %d" %(myname,i)
            
        s1,s2,s3,s4 = read_spectra_dada(firstSpectra, ifp,s1,s2,s3,s4)


        # if foff is negative, write the frequency channels in reverse order
        if foff < 0:   
            s1 = s1[::-1] # flip 1d array
            s2 = s2[::-1] # flip 1d array
            s3 = s3[::-1] # flip 1d array
            s4 = s4[::-1] # flip 1d array

        write_data(ofp1, s1)
        write_data(ofp2, s2)
        write_data(ofp3, s3)
        write_data(ofp4, s4)

        firstSpectra = 0
        
    ifp.close()
    ofp1.close()
    ofp2.close()
    ofp3.close()
    ofp4.close()

    t1 = time.time()
    total = t1-t0
    print "%s: time = %f", (myname,total)

    return



def send_string(ofp, str):
    ofp.write(struct.pack("<i",len(str)))
    ofp.write(str)

def send_int(ofp, str, value):
    send_string(ofp, str)
    ofp.write(struct.pack("<i",value))

def send_double(ofp, str, value):
    send_string(ofp, str)
    ofp.write(struct.pack("<d",value))
    

def open_file(fn):
    myname = "open_file"
    nch = proc_params.NCH
    if (proc_params.DEBUG):
        print("%s: nch = %d" %(myname, nch))

    # allocate mem for 4 spectra only
    s1 = np.zeros(nch, dtype=np.double);
    s2 = np.zeros(nch, dtype=np.double);
    s3 = np.zeros(nch, dtype=np.double);
    s4 = np.zeros(nch, dtype=np.double);

    fp = open(fn,'rb');
    return fp, s1, s2, s3, s4

def calc_mjd1(datestr):
    """
    Note: This is integer division ON PURPOSE (for once)
    """
    # this one just returns the integer part of mjd
    # example: 150720 
    year  = int(datestr[0:2]) + 2000
    month = int(datestr[2:4])
    day   = int(datestr[4:6])

    a = (14 - month)/12
    y = year + 4800 - a
    m = month + 12 * a - 3

    jdn = day + (153 * m + 2)/5 + 365 * y +y/4 - y/100 + y/400 - 32045
    mjd = jdn - 2400000.5

    return int(mjd)

def calc_mjd2(datestr):
    # example: 150720 17 40 00
    year  = int(datestr[0:2]) + 2000
    month = int(datestr[2:4])
    day   = int(datestr[4:6])

    a = (14 - month)/12
    y = year + 4800 - a
    m = month + 12 * a - 3

    jdn = day + (153 * m + 2)/5 + 365 * y +y/4 - y/100 + y/400 - 32045
    mjd = jdn - 2400000.5

    hh=datestr[7:9]; mm=datestr[10:12]; ss=datestr[13:15]
    hh=float(hh); mm=float(mm); ss=float(ss);
    mjd_frac = (hh + (mm/60.0) + (ss/3600.0))/24.0
    mjd = int(mjd) + mjd_frac

    return mjd


def get_nspec_dada(fn):

    myname = "get_nspec_dada"
    if (proc_params.DEBUG):
        print("%s: fn: %s" %(myname, fn))

    # number of bytes in one spectra
    # nch * 4 IF bands per file * 2 bytes per value (16 bit)
    ssize = (proc_params.NCH * 4 * 2) 

    statinfo = os.stat(fn)
    fsize = statinfo.st_size

    fsize = fsize - (512*8)  # account for 512 8-byte header

    nspectra = fsize / ssize

    if (proc_params.DEBUG):
        print "get_nspec_dada: nspectra = ", nspectra

    return nspectra


def write_header(ofn, fch1, source_name, tsamp, mjd):
    myname = "write_header"
    if (proc_params.DEBUG):
        print("%s: ofn: %s" %(myname, ofn))
        
    foff     = proc_params.FOFF
    nchans   = proc_params.NCH

    ofp = open(ofn, "wb")
    send_string(ofp, "HEADER_START")

    send_string(ofp, "rawdatafile")
    send_string(ofp, ofn)

    send_string(ofp, "source_name")
    send_string(ofp, source_name)

    send_int(ofp, "telescope_id",telescope_id)
    send_int(ofp, "machine_id",machine_id)
    send_int(ofp, "data_type",1) # filterbank
    send_double(ofp, "fch1",fch1)
    send_double(ofp, "foff",foff)
    send_int(ofp, "nchans",nchans)
    send_int(ofp, "nbits",obits)
    send_double (ofp, "tstart",mjd) 
    send_double(ofp, "tsamp",tsamp)
    send_int(ofp, "nifs",nifs)
    send_int(ofp, "barycentric",barycentric)

    send_string(ofp, "HEADER_END")

    ofp.close()


def write_data(ofp, data):
    fmt = '<' + str(len(data)) + 'f'
    obuf = struct.pack(fmt, *data)
    ofp.write(obuf)



def read_spectra_dada(first,fp, s1,s2,s3,s4):
    myname = "read_spectra_dada"
    nch = proc_params.NCH
    
    dat1_1 = np.zeros(nch/4, dtype=np.double);
    dat1_2 = np.zeros(nch/4, dtype=np.double);
    dat1_3 = np.zeros(nch/4, dtype=np.double);
    dat1_4 = np.zeros(nch/4, dtype=np.double);

    dat2_1 = np.zeros(nch/4, dtype=np.double);
    dat2_2 = np.zeros(nch/4, dtype=np.double);
    dat2_3 = np.zeros(nch/4, dtype=np.double);
    dat2_4 = np.zeros(nch/4, dtype=np.double);

    dat3_1 = np.zeros(nch/4, dtype=np.double);
    dat3_2 = np.zeros(nch/4, dtype=np.double);
    dat3_3 = np.zeros(nch/4, dtype=np.double);
    dat3_4 = np.zeros(nch/4, dtype=np.double);

    dat4_1 = np.zeros(nch/4, dtype=np.double);
    dat4_2 = np.zeros(nch/4, dtype=np.double);
    dat4_3 = np.zeros(nch/4, dtype=np.double);
    dat4_4 = np.zeros(nch/4, dtype=np.double);

    if (first):
        counter  = np.fromfile(fp, dtype=np.dtype('>Q'), count=512)

    dat1 = np.fromfile(fp, dtype=np.dtype('>u2'), count = nch)
    dat2 = np.fromfile(fp, dtype=np.dtype('>u2'), count = nch)
    dat3 = np.fromfile(fp, dtype=np.dtype('>u2'), count = nch)
    dat4 = np.fromfile(fp, dtype=np.dtype('>u2'), count = nch)

    #print "discard = ", np.uint64(discard)
        
    j=0
    for i in range(nch/4):
        dat1_1[i] = dat1[j]
        dat1_2[i] = dat1[j+1]
        dat1_3[i] = dat1[j+2]
        dat1_4[i] = dat1[j+3]

        dat2_1[i] = dat2[j]
        dat2_2[i] = dat2[j+1]
        dat2_3[i] = dat2[j+2]
        dat2_4[i] = dat2[j+3]

        dat3_1[i] = dat3[j]
        dat3_2[i] = dat3[j+1]
        dat3_3[i] = dat3[j+2]
        dat3_4[i] = dat3[j+3]

        dat4_1[i] = dat4[j]
        dat4_2[i] = dat4[j+1]
        dat4_3[i] = dat4[j+2]
        dat4_4[i] = dat4[j+3]

        j=j+4

    s1 = [dat1_1, dat1_2, dat1_3, dat1_4]
    s2 = [dat2_1, dat2_2, dat2_3, dat2_4]
    s3 = [dat3_1, dat3_2, dat3_3, dat3_4]
    s4 = [dat4_1, dat4_2, dat4_3, dat4_4]

    s1 = np.reshape(s1, nch)
    s2 = np.reshape(s2, nch)
    s3 = np.reshape(s3, nch)
    s4 = np.reshape(s4, nch)

    return s1,s2,s3,s4


def parse_dada_name(dada_fn):

    dirname = os.path.dirname(dada_fn)
    basename = os.path.basename(dada_fn)

    year = basename[0:4]
    mon  = basename[5:7]
    day  = basename[8:10]
    hh   = basename[11:13]
    mm   = basename[14:16]
    ss   = basename[17:19]

    tmp1 = basename.split("_")
    tmp2 = tmp1[1].split(".")
    bytes = int(tmp2[0])

    nsamps = bytes/2/proc_params.NCH
    year_2k = str(int(year) - 2000)

    datestr = year_2k+mon+day+" "+hh+" "+mm+" "+ss
    mjd = calc_mjd2(datestr)
    tmp3 = basename.split("_")
    ofn_base = tmp3[0]
    
    return (mjd, nsamps,ofn_base)

def get_dada_mjd(dada_fn):

    dirname = os.path.dirname(dada_fn)
    basename = os.path.basename(dada_fn)

    year = basename[0:4]
    mon  = basename[5:7]
    day  = basename[8:10]
    hh   = basename[11:13]
    mm   = basename[14:16]
    ss   = basename[17:19]

    tmp1 = basename.split("_")
    tmp2 = tmp1[1].split(".")
    bytes = int(tmp2[0])

    nsamps = bytes/2/proc_params.NCH
    year_2k = str(int(year) - 2000)

    datestr = year_2k+mon+day+" "+hh+" "+mm+" "+ss
    mjd = calc_mjd2(datestr)
    tmp3 = basename.split("_")
    ofn_base = tmp3[0]
    
    return (mjd)

def main2(argv):
    opts, args = getopt.getopt(argv, "hi:o:n:s:t:x:g:",["tsamp=", "source=", "tag="])
    # default values
    nproc = 0
    nskip = 0
    tsamp = 0
    tag = "default_tag"
    source = "default"
    idir = "/tmp/var"
    odir = "/tmp/var"

    msg_example= "usage: python proc_4ch -i idir -o odir --source b0329+54 --tsamp 513.024 --tag tag -n 5 -s 1"

    for opt, arg in opts:
        if opt == "-h":
            print "usage : ..."
            print msg_example
        elif opt == '-i':
            idir = arg
        elif opt == '-o':
            odir = arg
        elif opt in ("-t", "--tsamp"):
            tsamp = arg
        elif opt in ("-x", "--source"):
            source = arg
        elif opt in ("-g", "--tag"):
            tag = arg
        elif opt == '-n':
            nproc = arg
        elif opt == '-s':
            nskip = arg

    nproc = int(nproc)
    nskip = int(nskip)
    tsamp = float(tsamp)

    print "idir odir = ", idir, odir
    print "tsamp = ", tsamp
    print "tag = ", tag
    print "source = ", source
    print "nproc nskip = ", nproc, nskip

    tsamp    = tsamp * 1.0e-6 # convert usec to sec
    ch1freq  = proc_params.CH1FREQ
    ch2freq  = proc_params.CH2FREQ
    ch3freq  = proc_params.CH3FREQ
    ch4freq  = proc_params.CH4FREQ
    ch1label = proc_params.CH1LABEL
    ch2label = proc_params.CH2LABEL
    ch3label = proc_params.CH3LABEL
    ch4label = proc_params.CH4LABEL
    bw       = proc_params.BW

    filePattern = idir + '/*.dada'
    files = sorted(glob.iglob(filePattern))
    first = 1

    nfiles_processed = 0
    nToSkip = nskip

    for fn in files:
        print 'processing %s ...' % (fn)
        if nToSkip > 0:
            nToSkip = nToSkip - 1
        else:
            if first:
                # write headers
                mjd, nsamps, bname = parse_dada_name(fn)
                print("main: mjd %f nsamps %d bname %s" %(mjd, nsamps, bname))

                bname = tag + "-" + bname

                # set the names of the filterbank files using the first dada filename
                ofn1 = odir + "/" + bname + "_" + ch1label + ".fil"
                ofn2 = odir + "/" + bname + "_" + ch2label + ".fil"
                ofn3 = odir + "/" + bname + "_" + ch3label + ".fil"
                ofn4 = odir + "/" + bname + "_" + ch4label + ".fil"

                # BW = 1000 # MHz ==> channel bw is then 976562.5 Hz or slightly less than 1 MHz
                write_header(ofn1, ch1freq+bw, source, tsamp, mjd)
                write_header(ofn2, ch2freq+bw, source, tsamp, mjd)
                write_header(ofn3, ch3freq+bw, source, tsamp, mjd)
                write_header(ofn4, ch4freq+bw, source, tsamp, mjd)

            # split into 4 files
            dada_to_4ch(fn, ofn1, ofn2, ofn3, ofn4)
            nfiles_processed = nfiles_processed + 1
            first = 0
            if (nproc != 0) and (nproc == nfiles_processed):  # process only NPROC files
                break

msg_usage = "usage: python proc_4ch idir odir source dt (usec) [nproc]"
msg_example= "example: python proc_4ch /data2/15a201/s05-b0833 /data2/tmp b0329+54 513.024 tag"
msg_example= "example: python proc_4ch /data2/15a201/s05-b0833 /data2/tmp b0329+54 513.024 tag 5"

if __name__ == '__main__':
    from sys import argv

    #if len(sys.argv) < 6: print msg_usage; print msg_example; sys.exit(0)
    #if len(sys.argv) == 6: 
    #    nproc = 0
    #elif len(sys.argv) == 7: 
    #    nproc = int(argv[6])
    #main(argv[1], argv[2], argv[3], float(argv[4]), argv[5].rstrip(), nproc)

    main2(argv[1:])

