# __Date Last Revision__
# : 05/24/2017
raw_folder = '/data1/davidsw/sample_data/raw/'
processed_folder = '/data1/davidsw/sample_data/processed/'
scan_table_fp = '/home/davidsw/scan.table.all'
scan_log_fp = '/home/davidsw/mdscc.mars'

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
CH2TELESCOPE = 'DSS-54' # Channel 2 telescope
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
