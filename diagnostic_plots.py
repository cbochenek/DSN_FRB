import numpy as np,struct,time,datetime,sys,socket,os
import matplotlib.pyplot as plt
from sigpyproc.Readers import FilReader
import argparse
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from collections import deque
import scipy.ndimage as snd

def process_spec(spec):
	spec_reverse = np.flip(spec,0)
	#spec_flip4 = np.ravel(np.fliplr(spec_reverse.reshape(2048/4,4)))
	return spec_reverse

def process_data_arr(data_arr):
	#data_arr_good = np.flip(np.flip(data_arr.T.reshape(nsamples,2048/4,4),2).reshape(nsamples,2048),1)
	data_arr_good = np.flip(data_arr,1)
	return data_arr_good

def openfile(filename):
        files = FilReader(filename)
        nsamples = files.header['nsamples']
        data_arr = files.readBlock(0,nsamples)
        chans = np.array(range(files.header['nchans']))
        bandwidth = files.header['fch1']-files.header['fbottom']
        freqs = files.header['fbottom'] + chans/float(files.header['nchans'])*bandwidth


def dedisperse(data_arr,DM,tsamp,fs):
	#fs = np.linspace(1.28,1.53,2048)
	fs = fs/1000.
	dts = (0.00415/tsamp*DM*(fs[0]**-2-(fs)**-2)).astype(int)
	out_data_arr = np.zeros(data_arr.shape)
	for i in range(len(dts)):
		x = np.roll(data_arr[i],dts[len(fs)-1-i])
		out_data_arr[i] = x
	return out_data_arr
parser = argparse.ArgumentParser(description='Tell me what to do.')
parser.add_argument('-s', action = 'store_true',help='Plot the average spectrum')
parser.add_argument('-t', action = 'store_true', help = 'Plot power vs time summed over frequency')
parser.add_argument('-st', action = 'append', help = 'Plot power vs time in a specified frequency bin')
parser.add_argument('-im', action = 'store_true', help = 'Plot the dynamic spectrum')
parser.add_argument('-imd', action = 'store_true', help = 'Plot the dedispersed dynamic spectrum')
parser.add_argument('-sub', action = 'store_true', help = 'Compare spectrum with first light spectrum')
parser.add_argument('-HI', action = 'store_true', help = 'Save PNGs of 10 minute integrations of data around the H1 line.')
parser.add_argument('-fft', action = 'append', help ="FFT this channel")
parser.add_argument('-f', action = 'append', help = "Candidate number to open")
parser.add_argument('-b', action = 'append', help = "Plot the baseline with this resolution in seconds")
args = parser.parse_args()


cand_idx = int(args.f[0])
hd_cands = np.genfromtxt('/home/cbochenek/heimdall.cand', dtype=str)
hd_cand = hd_cands[np.where(hd_cands[:,1].astype(int)==cand_idx)[0][0]]
print hd_cand
boxcar_filters = np.array([1,2,3,5,10])
SNR = float(hd_cand[0])
ind = int(hd_cand[1])
cand_time_s = float(hd_cand[2])
cand_sample_width = boxcar_filters[int(hd_cand[3])]
cand_DM_trial_n = int(hd_cand[4])
cand_DM = float(hd_cand[5])
cand_members = int(hd_cand[6])
cand_data = hd_cand[7]
files = FilReader("%s" %(cand_data))
nsamples = files.header['nsamples']
data_arr = files.readBlock(ind-400,800)
tsamp = files.header['tsamp']
chans = np.array(range(files.header['nchans']))
bandwidth = files.header['fch1']-files.header['fbottom']
freqs = files.header['fbottom'] + chans/float(files.header['nchans'])*bandwidth
#integrations = np.array(range(nsamples))
#plot the average spectrum if the -s flag is seen
if args.s:
	spec_raw = data_arr.mean(axis=1)
	spec = process_spec(spec_raw)
	#print spec.min()
	plt.figure(0)
	plt.plot(freqs[:-5],spec[5:], 'k-')
	#plt.ylim(10.,14.)
	plt.xlabel('Frequency [MHz]')
	plt.ylabel('Power')
if args.t:
	tseries = data_arr.sum(axis = 0)
	plt.figure(1)
	plt.plot(tseries,'k-')
	plt.xlabel('Integration')
	plt.ylabel('Power')
if args.st !=None:
	spec_bins = np.array(args.st).astype(int)
	tseries_spec_bins = np.array([])
	for i in range(n_10m_integrations):
		data_arr_small = files.readBlock(i*286102,286102)
		int_spectra_raw = data_arr_small.mean(axis=1)
		int_spectra = process_spec(int_spectra_raw)
		spec_bins_power = int_spectra[spec_bins].sum()
		tseries_spec_bins = np.append(tseries_spec_bins,spec_bins_power)
	#tseries_spec_bins = data_arr[spec_bins,:].sum(axis=0)
	plt.figure(2)
	plt.plot(tseries_spec_bins,'k-')
	title = 'Bins '
	for arg in args.st:
		title = title +  arg + ' '
	plt.title(title)
	plt.ylabel('Power')
	plt.xlabel('Integration')
if args.im:
	#data_arr_good = process_data_arr(data_arr)
	data_arr_good = data_arr
	fig = plt.figure(3)
	gs = gridspec.GridSpec(2,2, width_ratios=[3,1], height_ratios=[1,3])
	ax0=plt.subplot(gs[2])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[3])
	data_arr_conv = snd.uniform_filter1d(data_arr_good,cand_sample_width,axis=1)
	data_arr_snr_dedisp = dedisperse(data_arr_conv,cand_DM,tsamp,freqs)
	print cand_DM
	ax0.imshow(data_arr_conv, aspect = 'auto', cmap = 'hot')
	ax0.set_title("%s" %cand_data,fontsize=12)
	ax0.set_xlabel("Integration \n%s %s %s %s %s %s %s" %(SNR, cand_idx, cand_time_s, cand_sample_width, cand_DM_trial_n, cand_DM,cand_members),fontsize=12)
	ax0.set_ylabel("Frequency Channel")
	ax2.plot(data_arr_conv.mean(axis=1),np.array(range(len(freqs))),color='k')
	ax2.set_ylim(0,len(freqs))
	ax1.plot(np.array(range(800)),data_arr_conv.mean(axis=0),color='k')
	ax2.set_ylim(np.max(np.array(range(len(freqs)))),np.min(np.array(range(len(freqs)))))
	ax1.set_xlim(np.min(np.array(range(800))),np.max(np.array(range(800))))
	#plt.figure(56)
	#plt.plot(range(len(data_arr_snr_dedisp.mean(axis=0))),data_arr_snr_dedisp.mean(axis=0),'k-')
if args.imd:
        #data_arr_good = process_data_arr(data_arr)
        data_arr_good = data_arr
        fig = plt.figure(4)
        gs = gridspec.GridSpec(2,2, width_ratios=[3,1], height_ratios=[1,3])
        ax0=plt.subplot(gs[2])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[3])
        data_arr_conv = snd.uniform_filter1d(data_arr_good,cand_sample_width,axis=1)
        data_arr_snr_dedisp = dedisperse(data_arr_conv,cand_DM,tsamp,freqs)
        print cand_DM
        ax0.imshow(data_arr_snr_dedisp, aspect = 'auto', cmap = 'hot')
	ax0.set_title("%s" %cand_data,fontsize=12)
        ax0.set_xlabel("Integration \n%s %s %s %s %s %s %s" %(SNR, cand_idx, cand_time_s, cand_sample_width, cand_DM_trial_n, cand_DM,cand_members),fontsize=12)
        ax0.set_ylabel("Frequency Channel")
        ax2.plot(data_arr_snr_dedisp.mean(axis=1),np.array(range(len(freqs))),color='k')
        ax2.set_ylim(0,len(freqs))
        ax1.plot(np.array(range(800)),data_arr_snr_dedisp.mean(axis=0),color='k')
        ax2.set_ylim(np.max(np.array(range(len(freqs)))),np.min(np.array(range(len(freqs)))))
        ax1.set_xlim(np.min(np.array(range(800))),np.max(np.array(range(800))))
if args.b != None:
	baseline = snd.uniform_filter1d(data_arr.mean(axis=0),int(float(args.b[0])/tsamp))
	plt.figure(5)
	plt.plot(baseline)
	plt.title("%s" %cand_data,fontsize=12)
        plt.xlabel("Integration \n%s %s %s %s %s %s %s" %(SNR, cand_idx, cand_time_s, cand_sample_width, cand_DM_trial_n, cand_DM,cand_members),fontsize=12)
	plt.ylabel("Power")

if args.HI:
        for i in range(64,65):
                data_arr_small = files.readBlock(i*15259,15259)
		data_arr_bp = data_arr_small[:,:100].mean(axis=0)
		#bp = np.ones(2048,15259)
                #int_spectra_raw = data_arr_small.mean(axis=1)
                #int_spectra = process_spec(int_spectra_raw)	
		#plt.plot(freqs,int_spectra)
		#plt.ylim(0,2500)
		plt.imshow(data_arr_small, aspect = 'auto', cmap = 'hot')
		plt.xlabel("Frequency [MHz]")
		plt.ylabel("Power")
		#plt.savefig("/home/user/peryton_search_%s.png" %i)
		#plt.clf()
if args.fft != None:
	spec_bins = np.array(args.fft).astype(int)
	count=0
	for chan in spec_bins:
		tseries = data_arr[chan]
		chan_fft = np.fft.fft(tseries)
		freqs_fft = np.fft.fftfreq(len(tseries),tsamp)
		chan_power = chan_fft*np.conj(chan_fft)
		plt.figure(10+count)
		plt.plot(freqs_fft[:len(tseries)/2],chan_power[:len(tseries)/2],'k-')
                plt.title("Channel %s" %chan)
                plt.xlabel("Frequency [Hz]")
                plt.ylabel("Power")
                plt.xscale('log')
                plt.yscale('log')	

plt.show()
