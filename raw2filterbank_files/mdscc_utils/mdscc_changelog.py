# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 12:41:35 2018

@author: sarabsha
"""
import datetime as dt

date = dt.datetime.strptime(datatag, '%ym%j').date()  # analysis:ignore
change_date_datatag = '16m329'
change_date = dt.datetime.strptime(change_date_datatag, '%ym%j').date()
if date > change_date:
    CH1FREQ = 2000.0 + BW  # Channel 1 center frequenccy (MHz) analysis:ignore
    CH2FREQ = 2000.0 + BW  # Channel 2 center frequenccy (MHz) analysis:ignore
    CH3FREQ = 1325.0 + BW  # Channel 3 center frequenccy (MHz) analysis:ignore
    CH4FREQ = 8100.0 + BW  # Channel 4 center frequenccy (MHz) analysis:ignore
    CH1LABEL = "S-LCP (Narrow L)"   # Channel 1 band name
    CH2LABEL = "S-RCP (Narrow L)"   # Channel 2 band name
    CH3LABEL = "L-LCP"              # Channel 3 band name
    CH4LABEL = "X-LCP"              # Channel 4 band name
    CH1TELESCOPE = 'DSS-63'  # Channel 1 telescope
    CH2TELESCOPE = 'DSS-63'  # Channel 2 telescope
    CH3TELESCOPE = 'DSS-63'  # Channel 3 telescope
    CH4TELESCOPE = 'DSS-63'  # Channel 4 telescope

    '''
    # Channels 400 to 680:
    truncate_freqs['ch1'] = [2000.0 - (400.0*FOFF), 2000.0 - (680.0*FOFF)]  # analysis:ignore
    # Channels 365 to 900:
    truncate_freqs['ch2'] = [2000.0 - (365.0*FOFF), 2000.0 - (900.0*FOFF)]  # analysis:ignore
    # Channels 580 to 900:
    truncate_freqs['ch3'] = [1325.0 - (580.0*FOFF), 1325.0 - (900.0*FOFF)]  # analysis:ignore
    # Channels 100 to 1024:
    truncate_freqs['ch4'] = [8100.0 - (100.0*FOFF), 8100.0 - (1024.0*FOFF)]  # analysis:ignore
    '''