#!/usr/bin/env python2
"""
Created on Wed May 24 14:46:39 2017

:author: sarabsha
"""
##### Import Default Libraries
import numpy as np
import os, glob, time, datetime
import getopt, sys
from math import log10

# ##### Import User defined Libraries
from mdscc_raw2filterbank_funcs import *



if __name__=="__main__":
    usage = "Usage: mdscc_raw2filterbank -s <scan-table> -c <config-file> -t <data-tag> -i <input-path> -o <output-path> -n <numer-of-files>"
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:c:t:i:o:n:",
                                   ["help", "scan-table=", "config-file=", "data-tag=",
                                   "input-path=", "output-path=", "writeS1", "writeS2", "writeSK"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        print usage
        sys.exit(2)

    num_process = 'all'
    st = 'all'
    scan_outpath = ''
    config_fp = 'mdscc_rawparam.py'

    output_S1raw = 1
    output_S2raw = 0
    output_SKestimator = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print usage
            sys.exit()
        elif o in ("-s", "--scan-table"):
            scan_table_fp = a
        elif o in ("-c", "--config-file"):
            config_fp = os.path.normpath(a)
        elif o in ("-t", "--scan-tag"):
            st = a
        elif o in ("-i", "--input-path"):
            datatag = os.path.basename(os.path.normpath(a))
            inpath_base = os.path.dirname(os.path.normpath(a))
        elif o in ("-o", "--output-path"):
            scan_outpath = a
        elif o in ("-n"):
            num_process = a
        elif o in ("--writeS1"):
            output_S1raw = 1
        elif o in ("--writeS2"):
            output_S2raw = 1
        elif o in ("--writeSK"):
            output_SKraw = 1
        else:
            print o
            assert False, "unhandled option"

    T = 1

    output_flags = {}
    output_flags['S1raw']           = output_S1raw
    output_flags['S2raw']           = output_S2raw
    output_flags['SKestimator']     = output_SKestimator

    execfile(config_fp)

    np.seterr(all='ignore')
    start_time = time.time()

    inpath = os.path.join(inpath_base, datatag)

    col_names = 'name,duration_min,duration_sec,date,start_hour,start_min,start_sec,tag'
    scan_table = np.genfromtxt(scan_table_fp, autostrip=True,
                         dtype=None,
                         names = col_names)

    if st == 'all':
        dt = int(datetime.datetime.strptime(datatag, '%ym%j').strftime('%y%m%d'))
        scan = scan_table[np.where(scan_table['date']==dt)[0]]
    else:
        scan = scan_table[np.where(scan_table['tag']==int(st))[0]]

    for i in range(len(scan)):
        print 'Begin processing files.'
        source_info = scan[i]

        tstart = get_tsart(source_info)

        source_name = source_info['name']
        scan_tag = str(source_info['tag'])

        infiles = sorted(glob.glob(os.path.join(inpath,'*'+(2-int(log10(float(scan_tag))))*'0'+scan_tag)))
        if infiles == []:
            raise IOError("No files to read!")
        bname = os.path.basename(infiles[0]).split("_")[0]

        print "scan date = ", datetime.datetime.strptime(str(source_info['date']),'%y%m%d').strftime('%m-%d-%y')
        print "source name = ", source_name
        print "scan tag = ", scan_tag
        print "bname = ", bname

        header['source_name'] = source_name
        header['tstart'] = tstart
        header['scan_tag'] = scan_tag
        header['bname'] = bname


        if st == 'all':
            scan_outpath = os.path.join(scan_outpath,scan_tag.zfill(2)+'_'+source_info['name'])
            if not os.path.exists(scan_outpath):
                os.makedirs(scan_outpath)
        else:
            if not os.path.exists(scan_outpath):
                sys.exit("Output directory doesn't exist!")


        M = float(NACCUM)
        if len(infiles)< T:
            T = len(infiles)

        f_scale = process_SK_ScaleCalc(infiles[T-1], header, M)

        if num_process != 'all':
            infiles = infiles[:int(num_process)]

        process_raw2filterbank(infiles, scan_outpath, header, f_scale, FCH1, M, output_flags)

        scan_outpath = ''
        print 'End processing ' + source_name+ ' files.'
        print "process time = %s seconds" %(time.time()-start_time)

    print 'End processing ' + datatag+ ' files.'
    print "process time = %s seconds" %(time.time()-start_time)