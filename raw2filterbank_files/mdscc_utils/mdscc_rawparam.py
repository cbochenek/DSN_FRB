# __Date Last Revision__
# : 05/22/2018

#################################
#
# PATHS:
#
#################################
raw_folder = '/data1/davidsw/sample_data/raw/'
processed_folder = '/data1/cbochenek/'
scan_table_fp = '/home/cbochenek/raw2filterbank_files/scan.table.all'
scan_log_fp = '/home/cbochenek/raw2filterbank_files/mdscc.mars'

#################################
#
# FPGA parameters:
#
#################################
NCH      = 1024         # number of spectral channels
SAMPRATE = 950.0        # Sampling rate (MS/s)
BW       = SAMPRATE/2.0 # Bandwidth (MHz)
FOFF     = -BW/NCH      # Width of the frequency channel
NACCUM   = 1024         # Number of accumulaitons
TSAMP    = 1e-6*NACCUM*\
            NCH/BW      # Sampling time (us)
# 1 accumulation = 1.024e-6 s = 1 * 1024 / 1e9
CH1FREQ  = 8100.0+BW    # Channel 1 center frequenccy (MHz)
CH2FREQ  = 8100.0+BW    # Channel 2 center frequenccy (MHz)
CH3FREQ  = 8100.0+BW    # Channel 3 center frequenccy (MHz)
CH4FREQ  = 8100.0+BW    # Channel 4 center frequenccy (MHz)
CH1LABEL = "X-LCP"# "S-LCP"# Channel 1 band name
CH2LABEL = "X-LCP"# "S-RCP"# Channel 2 band name
CH3LABEL = "X-RCP"# "L-LCP"# Channel 3 band name
CH4LABEL = "X-LCP"# "X-LCP"# Channel 4 band name
CH1TELESCOPE = 'DSS-55' # Channel 1 telescope
CH2TELESCOPE = 'DSS-65' # Channel 2 telescope
CH3TELESCOPE = 'DSS-63' # Channel 3 telescope
CH4TELESCOPE = 'DSS-63' # Channel 4 telescope

NSKIP    = 0
DEBUG    = 0            # debug flag


# set default header info here
telescope_id = 63
machine_id = 999

data_type = 1           # filterbank
obits = 32              # Output file bits

nifs = 1
barycentric = 0

#################################

mars_sats = ['M01O', 'MRO', 'MVN', 'MSL', 'MEX']

#################################
header = {}
header['telescope_id'] = telescope_id
header['machine_id'] = machine_id
header['data_type'] = data_type
header['foff'] = FOFF
header['nchans'] = NCH
header['obits'] = obits
header['tsamp'] = TSAMP
header['nifs'] = nifs
header['barycentric'] = barycentric

CHLABELS = {}
CHLABELS['ch1'] = CH1LABEL
CHLABELS['ch2'] = CH2LABEL
CHLABELS['ch3'] = CH3LABEL
CHLABELS['ch4'] = CH4LABEL

FCH1 = {}
FCH1['ch1'] = CH1FREQ
FCH1['ch2'] = CH2FREQ
FCH1['ch3'] = CH3FREQ
FCH1['ch4'] = CH4FREQ


ch_config = {}
ch_config['ch1'] = CH1TELESCOPE
ch_config['ch2'] = CH2TELESCOPE
ch_config['ch3'] = CH3TELESCOPE
ch_config['ch4'] = CH4TELESCOPE


mask_satfreqs = {}
mask_satfreqs['ch1'] = [8400.0, 8465.0]  # 8450]
mask_satfreqs['ch2'] = [8400.0, 8465.0]  # 8450]
mask_satfreqs['ch3'] = [8400.0, 8465.0]  # 8450]
mask_satfreqs['ch4'] = [8400.0, 8465.0]  # 8450]


truncate_freqs = {}
# Channels 680 to 811 (old: to 810 --> 8475.732422):
truncate_freqs['ch1'] = [8390.0, 8500.0]
#                       # [8380.17578125, 8506.81152344]
#                       # New frequencies 08/23/2017
#                       # [8415.429712, 8476.196289]  # [2660.5468, 2790.42956]
# Channels 101 to 1024 (old: from 100 --> 8146.386892):
truncate_freqs['ch2'] = [8146.85075900, 8573.14453197]
#                       # [2644.311455, 2892.4803]#
# Channels 101 to 1024 (old: from 100 --> 8146.386892):
truncate_freqs['ch3'] = [8146.85075900, 8573.14453197]
#                       # [2069.04286, 2217.4803]#
# Channels 101 to 1024 (old: from 100 --> 8146.386892):
truncate_freqs['ch4'] = [8146.85075900, 8573.14453197]  # 8575.000000

#################################
#
# Pulse search parameters:
#
#################################

search_threshold = 8.0
tw = 1.0 # (s) : extraction time window
