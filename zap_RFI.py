import numpy as np
import sigpyproc as sp
from sigpyproc.Readers import FilReader
from scipy import stats
import time
import scipy.ndimage as snd
import argparse
from scipy.signal import butter, lfilter
from scipy.signal import freqz,square
from scipy import fftpack

start_time = time.time()

parser = argparse.ArgumentParser(description='Tell me what to do.')
parser.add_argument('-f', action = 'append',help='Select the epoch to process')
parser.add_argument('-p', action = 'append',help='Choose the channels that contain 70m data')
parser.add_argument('-s', action = 'append',help='Select the scan number')
parser.add_argument('-d', action = 'append',help='Choose the base directory (e.g. /data3/cbochenek)')
args = parser.parse_args()

files = []
for arg in args.p:
	files.append(FilReader("%s/%s/mdscc_ch%s_S1_scan%s_60hzfilter.fil" %(args.d[0],args.f[0],arg,args.s[0])))
nsamples = files[0].header['nsamples']
GULP = int(250000)
NCHANS = files[0].header['nchans']
tsamp = files[0].header['tsamp']
BF = 6
BT = 40
bp_time = 1.  #seconds
bp_bins = int(bp_time/tsamp)
#bp_bins=0
smooth_time = .066  #seconds
smooth_bins = int(smooth_time/tsamp)
#smooth_bins = 0
old_bp = np.zeros(NCHANS)
bp_thresh = 2.5#/np.sqrt(bp_bins)
clipthresh_ss = clipthresh_ss = 64+16.
clipthresh_spec = clipthresh_spec = 64+2./np.sqrt(bp_bins)
clipthresh_ts = clipthresh_ts = 64.+40./np.sqrt(NCHANS)
clipthresh_block = (64)+16./np.sqrt(BF*BT)
zap = np.append(np.array(range(645,790)),np.array(range(985,1024))) 
old_bp = np.zeros(NCHANS)

def butter_bandpass_filter(data, lowcut, highcut, fs, order=6):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, data)
        return y
def butter_bandpass(lowcut, highcut, fs, order=6):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='bandstop')
        return b, a

def initialize_data(i, GULP, bp_bins,nsamples,files):
    if i==0:
        nsamp = GULP+bp_bins/2
        start = i*GULP
    elif i == nsamples/GULP and i != 0:
        nsamps_tail = nsamples-i*GULP
        nsamp = nsamps_tail+bp_bins/2
        start = i*GULP-bp_bins/2
    else:
        nsamp = GULP+bp_bins
        start = i*GULP-bp_bins/2
    data_arr = np.zeros((NCHANS,nsamp))
    norm = np.zeros((NCHANS,nsamp))
    zapped = np.zeros((NCHANS,nsamp))
    for j in range(len(args.p)):
        data_arr += files[j].readBlock(start,nsamp)
    zero_inds  = np.where(data_arr==0)
    data_arr[zero_inds[0],zero_inds[1]] = np.mean(data_arr)
    return nsamp, start, data_arr, norm, zapped

def correct_bp_and_zap_channels(data_arr, bp_bins, i, clipthresh_spec,bp_thresh,old_bp,norm, zapped):
        bp = (snd.uniform_filter1d(data_arr,bp_bins,axis=1))[:,bp_bins/2::bp_bins]
        for j in range(bp.shape[1]-1):
                div_bp = np.copy(bp[:,j])
                #div_bp[np.where(bp[:,j]==0)[0]] = 1e-5
                norm_bp = (((data_arr[:,j*bp_bins:(j+1)*bp_bins]).T/div_bp).T)*64.
                zapped_bp = zapped[:,j*bp_bins:(j+1)*bp_bins]
                if i == 0 and j == 0:
                        changed = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
                        norm_bp[changed] = 64.
                        zapped[changed,j*bp_bins:(j+1)*bp_bins] = 1.
                elif i!=0 and j == 0:
                        div_bp = np.copy(old_bp)
                        changed_1 = np.where(np.abs(bp[:,j]-old_bp)/div_bp > bp_thresh)[0]
                        changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
                        #div_bp[np.where(old_bp==0)[0]] = 1e-5
                        norm_bp[changed_1] = 64.
                        norm_bp[changed_2] = 64.
                        zapped[changed_1,j*bp_bins:(j+1)*bp_bins] = 1.
                        zapped[changed_2,j*bp_bins:(j+1)*bp_bins] = 1.
                else:
                        div_bp = np.copy(bp[:,j-1])
                        #div_bp[np.where(bp[:,j-1]==0)[0]] = 1e-5
                        changed_1 = np.where(np.abs(bp[:,j]-bp[:,j-1])/div_bp > bp_thresh)[0]
                        changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
                        norm_bp[changed_1] = 64.
                        norm_bp[changed_2] = 64.
                        zapped[changed_1,j*bp_bins:(j+1)*bp_bins] = 1.
                        zapped[changed_2,j*bp_bins:(j+1)*bp_bins] = 1.
                norm[:,j*bp_bins:(j+1)*bp_bins] = norm_bp
        div_bp1 = np.copy(bp[:,j+1])
        #div_bp1[np.where(bp[:,j+1]==0)[0]] = 1e-5
        div_bp0 = np.copy(bp[:,j])
        #div_bp0[np.where(bp[:,j]==0)[0]] = 1e-5
        norm_bp = (((data_arr[:,(j+1)*bp_bins:]).T/div_bp1).T)*64.
        changed_1 = np.where(np.abs(bp[:,j+1]-bp[:,j])/div_bp0 > bp_thresh)[0]
        changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
        norm_bp[changed_1] = 64.
        norm_bp[changed_2] = 64.
        zapped[changed_1,(j+1)*bp_bins:] = 1.
        zapped[changed_2,(j+1)*bp_bins:] = 1.
        norm[:,(bp.shape[1]-1)*bp_bins:] = norm_bp
        norm[zap] = 64.
        zapped[zap] = 1.
        #norm[227] = 64.
        #norm[52]=64.
        old_bp = np.copy(bp[:,-1])
        return old_bp, zapped, norm

def clean_pixels(norm_clean1, clipthresh_ss, zapped):
        clean1_inds = np.where(np.abs(norm_clean1 -64.)> clipthresh_ss-64.)
        norm_clean1[clean1_inds[0],clean1_inds[1]] = 64.
        zapped[clean1_inds[0],clean1_inds[1]] = 1.
        return zapped, norm_clean1

def clean_blocks(norm_clean2,BT,BF,NCHANS,nsamp,clipthresh_blocki,zapped):
        try:
                BT_pad = BT-norm_clean2.shape[1]%BT
                BF_pad = BF-norm_clean2.shape[0]%BF
                BT_BF_norm = np.pad(norm_clean2,((0,BF-norm_clean2.shape[0]%BF),(0,BT-norm_clean2.shape[1]%BT)),mode='constant',constant_values=64.).T.reshape(nsamp+BT_pad,(NCHANS+BF_pad)/BF,BF).mean(axis=2).T.reshape((NCHANS+BF_pad)/BF,(nsamp+BT_pad)/BT,BT).mean(axis=2)
                block_inds = np.where(np.abs(BT_BF_norm-64.) > clipthresh_block-64.)
                mask_arr = np.zeros(norm_clean2.shape)
                if block_inds[0].shape[0] != 0:
                        for j in range(block_inds[0].shape[0]):
                                mask_arr[block_inds[0][j]*BF:block_inds[0][j]*BF+BF,block_inds[1][j]*BT:block_inds[1][j]*BT+BT] = 1.
                zap_BF_BT = np.where(mask_arr == 1)
                norm_clean2[zap_BF_BT[0],zap_BF_BT[1]] = 64.
                zapped[zap_BF_BT[0],zap_BF_BT[1]] = 1.
        except ValueError:
                print("WARNING: BF and BT not set correctly. Skipping block filtering.")
        return zapped, norm_clean2

def clean_times(norm_clean3,zapped,clipthresh_ts,i):
        norm_tscrnch = norm_clean3.mean(axis=0)
        zap_ts = np.where(norm_tscrnch > clipthresh_ts)[0]
        norm_clean3[:,zap_ts] = 64.
        zapped[:,zap_ts] = 1.
        amount_zapped = float(zapped.sum())/(zapped.shape[0]*zapped.shape[1])
        print("Zapped %s percent of gulp %i." %(amount_zapped*100,i))
        return zapped, norm_clean3

def remove_baseline(norm_clean3,smooth_time,tsamp):
        baseline = snd.uniform_filter1d(norm_clean3.mean(axis=0),int(smooth_time/tsamp))
        baseline_subbed = norm_clean3 - baseline
        return baseline_subbed

def remove_10hz(baseline_subbed,zapped,tsamp):
        t = np.arange(baseline_subbed.shape[1]+46)
        boxcar = square(t/(1./2./np.pi/10./tsamp))
        av_pwr = baseline_subbed.mean(axis=0)
        A = np.fft.fft(boxcar[:-46])
        B = np.fft.fft(av_pwr)
        Br = B.conjugate()
        fftfreq = np.fft.fftfreq(B.shape[0],tsamp)
        tenHz_ind = np.where(np.abs(fftfreq-10.) == np.abs(fftfreq-10.).min())[0][0]
        if (B[tenHz_ind-5:tenHz_ind+5]*B[tenHz_ind-5:tenHz_ind+5].conjugate()).max() > 1e4:
                c = np.fft.ifft(A*Br)
                shift = np.argmax(c[:46])+1
                boxcar = boxcar[shift:shift+baseline_subbed.shape[1]]
                amp = np.mean(np.array([np.abs(np.median(av_pwr[np.where(boxcar>0)])),np.abs(np.median(av_pwr[np.where(boxcar<0)]))]))
                boxcar = amp*boxcar
                boxcar_arr = np.ones(baseline_subbed.shape)*boxcar
                dont_subtract = np.where(zapped==1)
                boxcar_arr[dont_subtract[0],dont_subtract[1]] = 0.
                no_boxcar = baseline_subbed - boxcar_arr
                baseline_no_boxcar = snd.uniform_filter1d(no_boxcar.mean(axis=0),int(smooth_time/tsamp))
                data_out = no_boxcar - baseline_no_boxcar
        else:
                data_out = baseline_subbed
        return data_out

def write_data(i,data_out,bp_bins,nsamples,GULP,args):
        if i == 0:
                out_data = data_out[:,:-bp_bins/2+1]
        elif i+1 == nsamples/GULP:
                out_data = data_out[:,bp_bins/2:]
        else:
                out_data = data_out[:,bp_bins/2:-bp_bins/2]
        outfile = sp.Filterbank.FilterbankBlock(out_data,files[0].header)
        outfile.toFile("%s/%s/mdscc_70m_S1_scan%s_60hzfilter_cleaned_%s.fil" %(args.d[0],args.f[0], args.s[0],i))
        return None

for i in range(nsamples/GULP+1):
	print(i)
	nsamp, start, data_arr, norm, zapped = initialize_data(i, GULP, bp_bins,nsamples,files)
	old_bp, zapped, norm = correct_bp_and_zap_channels(data_arr, bp_bins, i, clipthresh_spec,bp_thresh,old_bp,norm, zapped)
	zapped, norm_clean1 = clean_pixels(norm, clipthresh_ss, zapped)
	zapped, norm_clean2 = clean_blocks(norm_clean1,BT,BF,NCHANS,nsamp,clipthresh_block,zapped)
	zapped, norm_clean3 = clean_times(norm_clean2,zapped,clipthresh_ts,i)
	baseline_subbed = remove_baseline(norm_clean3,smooth_time,tsamp)
	print("Start remove 10hz")
	data_out = remove_10hz(baseline_subbed,zapped,tsamp)
	print("Done with 10hz nonsense")
	write_data(i,data_out,bp_bins,nsamples,GULP,args)
	"""if i==0:
		nsamp = GULP+bp_bins/2
		start = i*GULP
	elif i == nsamples/GULP and i != 0:
		nsamps_tail = nsamples-i*GULP
		nsamp = nsamps_tail+bp_bins/2
		start = i*GULP-bp_bins/2
	else:
		nsamp = GULP+bp_bins
		start = i*GULP-bp_bins/2
	data_arr = np.zeros((NCHANS,nsamp))
	norm = np.zeros((NCHANS,nsamp))
	zapped = np.zeros((NCHANS,nsamp))
	for j in range(len(args.p)):
		data_arr += files[j].readBlock(start,nsamp)
	zero_inds  = np.where(data_arr==0)
	data_arr[zero_inds[0],zero_inds[1]] = np.mean(data_arr)"""
	"""bp = (snd.uniform_filter1d(data_arr,bp_bins,axis=1))[:,bp_bins/2::bp_bins]
	for j in range(bp.shape[1]-1):
		div_bp = np.copy(bp[:,j])
		#div_bp[np.where(bp[:,j]==0)[0]] = 1e-5
		norm_bp = (((data_arr[:,j*bp_bins:(j+1)*bp_bins]).T/div_bp).T)*64.
		zapped_bp = zapped[:,j*bp_bins:(j+1)*bp_bins]
		if i == 0 and j == 0:
			changed = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
			norm_bp[changed] = 64.
			zapped[changed,j*bp_bins:(j+1)*bp_bins] = 1.
		elif i!=0 and j == 0:
			div_bp = np.copy(old_bp)
			changed_1 = np.where(np.abs(bp[:,j]-old_bp)/div_bp > bp_thresh)[0]
			changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
			#div_bp[np.where(old_bp==0)[0]] = 1e-5
			norm_bp[changed_1] = 64.
                        norm_bp[changed_2] = 64.
			zapped[changed_1,j*bp_bins:(j+1)*bp_bins] = 1.
			zapped[changed_2,j*bp_bins:(j+1)*bp_bins] = 1.
		else:
			div_bp = np.copy(bp[:,j-1])
			#div_bp[np.where(bp[:,j-1]==0)[0]] = 1e-5
			changed_1 = np.where(np.abs(bp[:,j]-bp[:,j-1])/div_bp > bp_thresh)[0]
			changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
			norm_bp[changed_1] = 64.
                        norm_bp[changed_2] = 64.
			zapped[changed_1,j*bp_bins:(j+1)*bp_bins] = 1.
			zapped[changed_2,j*bp_bins:(j+1)*bp_bins] = 1.
		norm[:,j*bp_bins:(j+1)*bp_bins] = norm_bp
	div_bp1 = np.copy(bp[:,j+1])
	#div_bp1[np.where(bp[:,j+1]==0)[0]] = 1e-5
	div_bp0 = np.copy(bp[:,j])
        #div_bp0[np.where(bp[:,j]==0)[0]] = 1e-5
	norm_bp = (((data_arr[:,(j+1)*bp_bins:]).T/div_bp1).T)*64.
	changed_1 = np.where(np.abs(bp[:,j+1]-bp[:,j])/div_bp0 > bp_thresh)[0]
	changed_2 = np.where(np.abs(norm_bp.mean(axis=1)-64.) > clipthresh_spec-64.)[0]
	norm_bp[changed_1] = 64.
	norm_bp[changed_2] = 64.
	zapped[changed_1,(j+1)*bp_bins:] = 1.
	zapped[changed_2,(j+1)*bp_bins:] = 1.
	norm[:,(bp.shape[1]-1)*bp_bins:] = norm_bp
	norm[zap] = 64.
	zapped[zap] = 1.
	#norm[227] = 64.
	#norm[52]=64.
	old_bp = np.copy(bp[:,-1])"""
	"""norm_clean1 = norm
	clean1_inds = np.where(np.abs(norm_clean1 -64.)> clipthresh_ss-64.)
	norm_clean1[clean1_inds[0],clean1_inds[1]] = 64.
	zapped[clean1_inds[0],clean1_inds[1]] = 1.
	norm_clean2 = norm_clean1"""
	#norm_fscrnch = norm_clean1.mean(axis=1)
	#norm_clean2[np.where(norm_fscrnch > clipthresh_spec)[0]] = 64.
	"""try:
		BT_pad = BT-norm_clean2.shape[1]%BT
		BF_pad = BF-norm_clean2.shape[0]%BF
		BT_BF_norm = np.pad(norm_clean2,((0,BF-norm_clean2.shape[0]%BF),(0,BT-norm_clean2.shape[1]%BT)),mode='constant',constant_values=64.).T.reshape(nsamp+BT_pad,(NCHANS+BF_pad)/BF,BF).mean(axis=2).T.reshape((NCHANS+BF_pad)/BF,(nsamp+BT_pad)/BT,BT).mean(axis=2)
		block_inds = np.where(np.abs(BT_BF_norm-64.) > clipthresh_block-64.)
		mask_arr = np.zeros(norm_clean2.shape)
		if block_inds[0].shape[0] != 0:
			for j in range(block_inds[0].shape[0]):
    				mask_arr[block_inds[0][j]*BF:block_inds[0][j]*BF+BF,block_inds[1][j]*BT:block_inds[1][j]*BT+BT] = 1.
		zap_BF_BT = np.where(mask_arr == 1)
		norm_clean2[zap_BF_BT[0],zap_BF_BT[1]] = 64.
		zapped[zap_BF_BT[0],zap_BF_BT[1]] = 1.
	except ValueError:
		print "WARNING: BF and BT not set correctly. Skipping block filtering."
	norm_clean3 = norm_clean2"""
	"""norm_tscrnch = norm_clean3.mean(axis=0)
	zap_ts = np.where(norm_tscrnch > clipthresh_ts)[0]
        norm_clean3[:,zap_ts] = 64.
	zapped[:,zap_ts] = 1.
	amount_zapped = float(zapped.sum())/(zapped.shape[0]*zapped.shape[1])
	print "Zapped %s percent of gulp %i." %(amount_zapped*100,i)"""
	#norm_clean4 = np.copy(norm_clean3)
	#kurtosis = stats.kurtosis(norm_clean4,axis=1)/(np.std(norm_clean4,axis=1)**2)
	#norm_clean4[np.where(kurtosis > kurt_cut)[0]] = 64.
	"""baseline = snd.uniform_filter1d(norm_clean3.mean(axis=0),int(smooth_time/tsamp))
	baseline_subbed = norm_clean3 - baseline"""
	"""t = np.arange(baseline_subbed.shape[1]+46)
        boxcar = square(t/(1./2./np.pi/10./tsamp))
        av_pwr = baseline_subbed.mean(axis=0)
        A = fftpack.fft(boxcar[:-46])
        B = fftpack.fft(av_pwr)
        Br = B.conjugate()
        fftfreq = np.fft.fftfreq(B.shape[0],tsamp)
        tenHz_ind = np.where(np.abs(fftfreq-10.) == np.abs(fftfreq-10.).min())[0][0]
        if (B[tenHz_ind-5:tenHz_ind+5]*B[tenHz_ind-5:tenHz_ind+5].conjugate()).max() > 1e4:
                c = fftpack.ifft(A*Br)
                shift = np.argmax(c[:46])+1
                boxcar = boxcar[shift:shift+baseline_subbed.shape[1]]
                amp = np.mean(np.array([np.abs(np.median(av_pwr[np.where(boxcar>0)])),np.abs(np.median(av_pwr[np.where(boxcar<0)]))]))
                boxcar = amp*boxcar
		boxcar_arr = np.ones(baseline_subbed.shape)*boxcar
		dont_subtract = np.where(zapped==1)
		boxcar_arr[dont_subtract[0],dont_subtract[1]] = 0.
                no_boxcar = baseline_subbed - boxcar_arr
                baseline_no_boxcar = snd.uniform_filter1d(no_boxcar.mean(axis=0),int(smooth_time/tsamp))
                data_out = no_boxcar - baseline_no_boxcar
        else:
                data_out = baseline_subbed"""
	"""if i == 0:
		out_data = data_out[:,:-bp_bins/2+1]
	elif i+1 == nsamples/GULP:
		out_data = data_out[:,bp_bins/2:]
	else:
		out_data = data_out[:,bp_bins/2:-bp_bins/2]
	outfile = sp.Filterbank.FilterbankBlock(out_data,files[0].header)
	outfile.toFile("/data3/cbochenek/%s/mdscc_70m_S1_scan%s_60hzfilter_cleaned_%s.fil" %(args.f[0], args.s[0],i))"""
print("--- %s seconds ---" %(time.time()-start_time))
