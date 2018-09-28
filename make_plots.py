import numpy as np,struct,time,datetime,sys,socket,os
import matplotlib.pyplot as plt
from sigpyproc.Readers import FilReader
import argparse
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
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


hd_cands = np.genfromtxt('/home/cbochenek/heimdall.cand', dtype=str)
boxcar_filters = np.array([1,2,3,5,10])
SNR = hd_cands[:,0].astype(float)
inds = hd_cands[:,1].astype(int)
cand_time = hd_cands[:,2].astype(float)
cand_sample_width = boxcar_filters[hd_cands[:,3].astype(int)]
cand_DM_trial_n = hd_cands[:,4].astype(int)
cand_DM = hd_cands[:,5].astype(float)
cand_members = hd_cands[:,6].astype(int)
cand_data = hd_cands[:,7]
good_inds = np.where(SNR >= 8.7)[0]

with PdfPages("candidates.pdf") as pdf:
	for ind in good_inds:
		files = FilReader("%s" %(cand_data[ind]))
		nsamples = files.header['nsamples']
		data_arr = files.readBlock(inds[ind]-150,300)
		tsamp = files.header['tsamp']
		chans = np.array(range(files.header['nchans']))
		bandwidth = files.header['fch1']-files.header['fbottom']
		freqs = files.header['fbottom'] + chans/float(files.header['nchans'])*bandwidth
		smooth_bins = int(0.066/tsamp)
		plt.figure(ind,figsize=(20,10))
		gs = gridspec.GridSpec(nrows = 6,ncols = 1)
		ax0=plt.subplot(gs[0])
		ax1=plt.subplot(gs[1])
		ax2=plt.subplot(gs[2])
		ax3=plt.subplot(gs[3])
		ax4=plt.subplot(gs[4])
		ax5=plt.subplot(gs[5])
		data_arr_good = data_arr[:,100:200]
		data_arr_conv = snd.uniform_filter1d(data_arr_good,cand_sample_width[ind],axis=1)
		data_arr_snr_dedisp = dedisperse(data_arr_conv,cand_DM[ind],tsamp,freqs)
		baseline = snd.uniform_filter1d(data_arr_good.mean(axis=0),smooth_bins)
		ax0.set_title("%s" %(cand_data[ind]),fontsize=12)
		ax5.set_xlabel("Integration \n%s %s %s %s %s %s %s" %(SNR[ind], inds[ind], cand_time[ind], cand_sample_width[ind], cand_DM_trial_n[ind], cand_DM[ind],cand_members[ind]),fontsize=12)
		ax0.plot(data_arr_conv.mean(axis=0))
		ax0.vlines([50,50+cand_sample_width[ind]],0,1,transform = ax0.get_xaxis_transform())
		ax2.plot(data_arr_snr_dedisp.mean(axis=0))
		ax2.vlines([50,50+cand_sample_width[ind]],0,1,transform = ax2.get_xaxis_transform())
		ax1.imshow(data_arr_conv,cmap='hot',aspect='auto')
		ax5.plot(np.array(range(len(freqs))),data_arr_conv.mean(axis=1))
		ax3.imshow(data_arr_snr_dedisp,cmap='hot',aspect='auto')
		ax4.plot(baseline)
		pdf.savefig()
		plt.close()
