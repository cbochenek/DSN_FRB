# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 10:44:12 2018

@author: davidsw
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import time
from threading import Thread
import timeit
import subprocess
import mdscc_rawparam
from subprocess import call, check_call, check_output, Popen
import filterbank
import filterbank_class as fb
import h5py
from scipy import signal
import multiprocessing
from multiprocessing import Queue
import math
import fb_merge as f1
# ----------------------------------------------------------
# ----------------------------------------------------------

#datatags = ['18m055']
scan_tag = 'all'
num_process = 'all'
tic = timeit.default_timer()
# mdscc_utils_path = os.path.join(os.path.dirname(sys.argv[0]), 'mdscc_utils')
mdscc_utils_path = '/home/cbochenek/raw2filterbank_files/' + \
                   'mdscc_utils/'
os.environ['PATH'] += os.pathsep + mdscc_utils_path
sys.path.append(mdscc_utils_path)

zaplist_fn = os.path.join(mdscc_utils_path, '60hz.10hz.zaplist')

# config_fp = '~/my_python_library/mdscc_dataprocess/mdscc_utils/' +  \
#             'mdscc_rawparam.py'
config_fp = os.path.join(mdscc_utils_path, 'mdscc_rawparam.py')
execfile(os.path.normpath(config_fp))
print config_fp
print "mdscc_utils_path %s" %mdscc_utils_path

flag = {'raw2filterbank':1,
        'lowres': 1,
        'truncate': 1,
        'mask': 0,
        'new_mask': 1,
        'rfifind': 1,
        'prepfil': 1,
        'prepdata': 1,
        'realfft': 1,
        'timeseries': 1,
        'bandpass': 1,
        'single_pulse_search': 0,
        }
if len(sys.argv) > 1:
    datatags = sys.argv[1:]
else:
    raise IOError("datatag is needed.")

for datatag in datatags:
    total_start = time.time()
    # ----------------------------------------------------------
    # Convert raw to filterbank
    # ----------------------------------------------------------
    input_p = os.path.join(raw_folder, datatag)
    output_p = os.path.join(processed_folder, datatag)
    in_p = input_p
    out_p = output_p

    convert_cmd = os.path.join(mdscc_utils_path, 'mdscc_raw2filterbank.py') \
        + ' -s ' + scan_table_fp \
        + ' -c ' + config_fp \
        + ' -t ' + scan_tag \
        + ' -i ' + input_p \
        + ' -o ' + output_p \
        + ' -n ' + num_process

    if flag['raw2filterbank'] == 1:
        print 'Converting raw files to filterbank files ...'
        print 'cmd>> ', convert_cmd
        sys.stdout.flush()
        start = time.time()
        subprocess.call(convert_cmd, shell=True)
        print "Done raw files to filterbank files in " + \
              "%s seconds" % (time.time() - start)
