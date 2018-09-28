import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.cm.jet.set_bad('darkgray',1.0)
#import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy import fftpack
import sys

sys.path.append('/home/davidsw/utils')
import filterbank_class as fb
#########################################
def pulse_search(fn, threshold):
    de_disp_fn = fn[:-4]+".dedisp"
    single_pulse_fn = de_disp_fn+".singlepulse"

    fn = fn + ".fil"

#    prepdata_cmd = "prepdata -filterbank  " + fn + " -o " + de_disp_fn + " -noclip -nobary -dm 0.0 "
#    print prepdata_cmd
#    os.system(prepdata_cmd)

    single_pulse_search_cmd = "single_pulse_search.py -b -m 100 -t " + str(threshold) + " " + de_disp_fn + ".dat"
    print single_pulse_search_cmd
    os.system(single_pulse_search_cmd)
    return single_pulse_fn

#########################################
def read_tw(fn, start_sample, nspectra, pulse_sample):
    if start_sample < 1:
        start_sample = 1
    ext_fn = fn[:-4]+'_'+str(pulse_sample)+".ext"
    extract_cmd = "extract " + fn + " " + str(start_sample) + " " + str(nspectra) + " > " + ext_fn
#    print extract_cmd
    os.system(extract_cmd)
    fbc = fb.Filterbank(ext_fn,1)
    data = fbc.rawData
    time =  np.arange(fbc.nsamps)*fbc.tsamp # (start_sample + np.arange(fbc.nsamps))*fbc.tsamp
    freq = np.arange(fbc.nchans)*fbc.foff + fbc.fch1
    return time, freq, data

#########################################
def setup_fig(labels, fig_title):
    fig = plt.figure(figsize=(12,12))
    fig.suptitle(fig_title, fontsize=14,y=0.97)
    gs = gridspec.GridSpec(2, 2)
    gs.update(left=0.1, right=0.925, wspace=0.35)
    axarr = []

    # Draw panels:
    for i in range(4):
        gss = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[i])
        params = gss.get_subplot_params().__dict__
        del params['validate']
        fig, axarr_tmp = plt.subplots(4, sharex=True, num=fig.number, squeeze = True,
                                    gridspec_kw=params)
        axarr.append(axarr_tmp)
        axarr[i][0].set_title(labels[i])
        axarr[i][0].set_ylabel('S1 ts')
        #axarr[i][0].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
        axarr[i][0].locator_params(axis='y', nbins=4)
        axarr[i][0].grid(linestyle='dashed', dash_joinstyle='round', dash_capstyle='round')

        axarr[i][1].set_ylabel('S1 freqs (GHz)')
        axarr[i][1].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
        axarr[i][1].locator_params(axis='y', nbins=4)
        axarr[i][1].grid(linestyle='dashed', dash_joinstyle='round', dash_capstyle='round')

        axarr[i][2].set_ylabel('SK ts')
        #axarr[i][2].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
        axarr[i][2].locator_params(axis='y', nbins=4)
        axarr[i][2].grid(linestyle='dashed', dash_joinstyle='round', dash_capstyle='round')

        axarr[i][3].set_ylabel('SK freqs (GHz)')
        axarr[i][3].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
        axarr[i][3].locator_params(axis='y', nbins=4)
        axarr[i][3].grid(linestyle='dashed', dash_joinstyle='round', dash_capstyle='round')
#        axarr[i][3].set_xlabel('time (s)')


    return fig, axarr
#########################################
class ScalarFormatterForceFormat(ticker.ScalarFormatter):
    def _set_format(self, other, another):
        self.format = "%1.1f"
sfmt = ScalarFormatterForceFormat(useLocale=True,useOffset=True)
sfmt.set_scientific(True)
sfmt.set_powerlimits((-2,3))
#########################################
def read_pulses(ref_fn, single_pulse_fn, filepath, tw, tsamp,labels, fig_title):
    nspectra = round(tw/tsamp)

#    de_disp_fn = ref_fn + ".dedisp"
#    single_pulse_fn = de_disp_fn + ".singlepulse"
    fig_fn = ref_fn

    ref_fn = ref_fn + ".fil"

    pulse_infos =  np.loadtxt(single_pulse_fn)
    if pulse_infos.size == 0:
        return
    if pulse_infos.ndim == 1:
        pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

    pn = 0
    for pulse_info in pulse_infos:
        pn += 1
#        if pn < 1946:
#            continue
#        if pn > 1946:
#            quit()
#        if pn > 1:
#            quit()
        pulse_time = pulse_info[2]
        pulse_sigma = pulse_info[1]
        pulse_width = pulse_info[4]
        if pulse_width > 0:
            print 'Begin producing pulse number ' + str(pn) +' figure.'

            start_sample = round((pulse_time-tw/2)/tsamp)

            pulse_fig_title = fig_title + '\n time = ' + str(pulse_time) + ' s ' + \
                ' , index = ' + str(int(pulse_info[3])) + \
                ' , sigma = ' + str(pulse_sigma) + \
                ' , width = ' + str(pulse_width)
            fig, axarr = setup_fig(labels, pulse_fig_title)


            chtags = ['ch1', 'ch2', 'ch3', 'ch4']
            for t in range(len(chtags)):
                ch = chtags[t]

                tag = "S1"
                fn = filepath + "_"  + ch + "_" + tag + "_new.fil"
                time, freq, data = read_tw(fn, start_sample, nspectra, int(pulse_info[3]))

                ts = data.sum(axis=1)
#                ts = ts - ts.min()
#                ts = np.ma.masked_less_equal(ts, 0)
#                ts = 10*np.log10(ts)
                axarr[t][0].plot(time, ts)

#                data = data - data.min()
#                data = np.ma.masked_less_equal(data, 0)
#                data = 10*np.log10(data)


                flat_data = data.flatten()
                vmin    = np.min(flat_data)
                vmax    = np.max(flat_data)
                vmean   = np.mean(flat_data)
                vmedian = np.median(flat_data)
                vstd    = np.std(flat_data)
                zmin = vmean - 3*vstd
                zmax = vmean + 3*vstd

                im = axarr[t][1].imshow(data.T, aspect = 'auto',
                                   extent=[min(time),max(time), min(freq)/1000.0, max(freq)/1000.0],
                                   vmin=zmin, vmax=zmax)


                rect = axarr[t][1].get_position()
                rect.x0 = rect.x1 + 0.005
                rect.x1 = rect.x0 + 0.015
                cax = plt.axes([rect.x0,rect.y0,rect.x1-rect.x0,rect.y1-rect.y0])
                cbar = fig.colorbar(im, cax=cax, ticks = [zmin, zmax], format=sfmt, extend='both')#, aspect=5, extend='both', pad=0.0, use_gridspec=False, anchor=(-10.0,0.0))
                cbar.ax.tick_params(labelsize=10)
                cbar.ax.yaxis.get_offset_text().set_fontsize(10)


#                tag = "SKts"
#                fn = filepath + "_"  + ch + "_" + tag + ".fil"
#                time, data = read_tw(fn, start_sample, nspectra, pn)
#                axarr[t][2].plot(time, data)

                tag = "SK"
                fn = filepath + "_"  + ch + "_" + tag + "_new.fil"
                time, freq, data = read_tw(fn, start_sample, nspectra, int(pulse_info[3]))

                ts = data.sum(axis=1)
#                ts = ts - ts.min()
#                ts = np.ma.masked_less_equal(ts, 0)
#                ts = 10*np.log10(ts)
                axarr[t][2].plot(time, ts)

#                data = data - data.min()
#                data = np.ma.masked_less_equal(data, 0)
#                data = 10*np.log10(data)


                flat_data = data.flatten()
#                if ch == 'ch1':
#                    flat_data = data[:,50:250].flatten()
                vmin    = np.min(flat_data)
                vmax    = np.max(flat_data)
                vmean   = np.mean(flat_data)
                vmedian = np.median(flat_data)
                vstd    = np.std(flat_data)
                zmin = vmean - 3*vstd
                zmax = vmean + 3*vstd
#                if ch == 'ch1':
#                    zmin = vmean - 1*vstd
#                    zmax = vmean + .8*vstd

                im = axarr[t][3].imshow(data.T, aspect = 'auto',
                                   extent=[min(time),max(time), min(freq)/1000.0, max(freq)/1000.0],
                                   #norm=colors.SymLogNorm(linthresh=1.5,vmin=zmin, vmax=zmax))
                                   vmin=zmin, vmax=zmax)

                rect = axarr[t][3].get_position()
                rect.x0 = rect.x1 + 0.005
                rect.x1 = rect.x0 + 0.015
                cax = plt.axes([rect.x0,rect.y0,rect.x1-rect.x0,rect.y1-rect.y0])
                cbar = fig.colorbar(im, cax=cax, ticks = [zmin, zmax], format=sfmt, extend='both')#, aspect=5, extend='both', pad=0.0, use_gridspec=False, anchor=(-10.0,0.0))
                cbar.ax.tick_params(labelsize=10)
                cbar.ax.yaxis.get_offset_text().set_fontsize(10)
                cbar.update_ticks()


                axarr[t][3].set_xlabel('Time (s)')

            fig.savefig(fig_fn+'_'+str(pn)+'.png')
            plt.close(fig)
            print 'End producing pulse number ' + str(pn) +' figure.'

            #break
            #raw_input("Press the <ENTER> key to continue...")
        #return

#########################################
def find_datatag_ref_ch(fp, mars_sats, ch_config):
    col_names = 'DATE,DAY,START,BOT,EOT,END,FACILITY,USER,ACTIVITY'.lower()

    scan_log = np.genfromtxt(fp, autostrip=True,
                             skip_header=6,
                             dtype=None, invalid_raise=False, #max_rows=140,
                             filling_values=-1,missing_values='',usecols=(0,1,2,3,4,5,6,7,8),
                             names = col_names)

#    days = np.unique(scan_log['day'])
    days = scan_log['day']
    refs = {}
    prev_day = 0
    year_i = 16
    for day in days:
        if day != -1:
            if day < prev_day:
                year_i += 1
            prev_day = day
            year = str(year_i)
#            if day < 329 :
#                year = '17'
#            else:
#                year = '16'
            f = np.unique(scan_log['facility'][np.where(scan_log['day'] == day)])[0]
            refs[year+'m'+str(day).zfill(3)] = [key for key, value in ch_config.iteritems() if f in value] #ch_config[f]
#    print sorted(refs.keys())

    return refs
#########################################
def pulse_scatterplot(fb_filepath, ref_chs, ref_tag):
    fig = plt.figure()
    plt.hold(True)
    plt.xlabel('time (s)')
    plt.ylabel('Channel #')
    plt.title('area: sigma, color: width')
    for ref_ch in ref_chs:

        ref_fn = fb_filepath + "_" + ref_ch + "_" + ref_tag #+ ".fil"

        de_disp_fn = ref_fn + ".dedisp"
        single_pulse_fn = de_disp_fn + ".singlepulse"
        fig_fn = fb_filepath+'_'+'scatterplot.png'

        ref_fn = ref_fn + ".fil"

        pulse_infos =  np.loadtxt(single_pulse_fn)
        if pulse_infos.size == 0:
            continue
        if pulse_infos.ndim == 1:
            pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

        pulse_times = pulse_infos[:,2]
        pulse_sigmas = pulse_infos[:,1]
        pulse_width = pulse_infos[:,4]
        area = pulse_sigmas #np.pi * ((pulse_sigmas - pulse_sigmas.min()))**2
        color = pulse_width

        pulse_ch = int(ref_ch[2])*np.ones(len(pulse_times))

        plt.scatter(pulse_times, pulse_ch, s=area, c=color, cmap = 'viridis', alpha = 0.5)

    if plt.gci():
        plt.colorbar()
    fig.savefig(fig_fn)
    plt.close(fig)
#########################################
def pulse_rateplot(fb_filepath, ref_chs, ref_tag,labels):
    fig = plt.figure()
    plt.hold(True)
    plt.xlabel('time (s)')
    plt.ylabel('number of pulses per bin')
    pulse_times = []#np.zeros()
    for ref_ch in ref_chs:

        ref_fn = fb_filepath + "_" + ref_ch + "_" + ref_tag #+ ".fil"

        de_disp_fn = ref_fn + ".dedisp"
        single_pulse_fn = de_disp_fn + ".singlepulse"
        fig_fn = fb_filepath+'_'+'pulserate.png'

        ref_fn = ref_fn + ".fil"

        pulse_infos =  np.loadtxt(single_pulse_fn)
        if pulse_infos.size == 0:
            continue
        if pulse_infos.ndim == 1:
            pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

        pulse_times.append(pulse_infos[:,2])

    plt.hist(pulse_times,25, label=[labels[int(ref_ch[2])-1] for ref_ch in ref_chs])#labels[int(ref_ch[2])-1])#, normed=1)#, facecolor='green')
    plt.legend(loc=2)

    fig.savefig(fig_fn)
    plt.close(fig)
#########################################
def pulse_SNRhistplot(fb_filepath, ref_chs, ref_tag,labels):
    fig = plt.figure()
    plt.hold(True)
    plt.xlabel('sigma')
    plt.ylabel('number of sigmas per bin')

    pulse_sigmas=[]
    for ref_ch in ref_chs:

        ref_fn = fb_filepath + "_" + ref_ch + "_" + ref_tag #+ ".fil"

        de_disp_fn = ref_fn + ".dedisp"
        single_pulse_fn = de_disp_fn + ".singlepulse"
        fig_fn = fb_filepath+'_'+'SNRhist.png'

        ref_fn = ref_fn + ".fil"

        pulse_infos =  np.loadtxt(single_pulse_fn)
        if pulse_infos.size == 0:
            continue
        if pulse_infos.ndim == 1:
            pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

        pulse_sigmas.append(pulse_infos[:,1])

    plt.hist(pulse_sigmas,25, label=[labels[int(ref_ch[2])-1] for ref_ch in ref_chs])
    plt.legend(loc=1)

    fig.savefig(fig_fn)
    plt.close(fig)
#########################################
def pulse_width_histplot(fb_filepath, ref_chs, ref_tag,labels):
    fig = plt.figure()
    plt.hold(True)
    plt.xlabel('width')
    plt.ylabel('number of width values per bin')

    pulse_width=[]
    for ref_ch in ref_chs:

        ref_fn = fb_filepath + "_" + ref_ch + "_" + ref_tag #+ ".fil"

        de_disp_fn = ref_fn + ".dedisp"
        single_pulse_fn = de_disp_fn + ".singlepulse"
        fig_fn = fb_filepath+'_'+'width_hist.png'

        ref_fn = ref_fn + ".fil"

        pulse_infos =  np.loadtxt(single_pulse_fn)
        if pulse_infos.size == 0:
            continue
        if pulse_infos.ndim == 1:
            pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

        pulse_width.append(pulse_infos[:,4])

    plt.hist(pulse_width,25, label=[labels[int(ref_ch[2])-1] for ref_ch in ref_chs])
    plt.legend(loc=1)

    fig.savefig(fig_fn)
    plt.close(fig)

#########################################
def spectrogram(fn):
    fbc = fb.Filterbank(fn,1)
    data = fbc.rawData - fbc.rawData.mean(axis=0)
    ts = fbc.tsamp
    fs = 1.0/(2.0*ts)
    Ns = fbc.nsamps
    Nf = int(Ns/2.0)
    Nch = fbc.nchans

    fftdata = abs(fftpack.fft(data,axis=0)[:Nf])**2
    fftfreq = fftpack.fftfreq(Ns, d=ts)[:Nf]
    freq = np.arange(fbc.nchans)*fbc.foff + fbc.fch1

    return fftfreq, freq, fftdata
#########################################
def plot_spectrograms(ref_fn, single_pulse_fn, filepath, tsamp, labels, fig_title):

    fig_fn = ref_fn

    ref_fn = ref_fn + ".fil"

    pulse_infos =  np.loadtxt(single_pulse_fn)
    if pulse_infos.size == 0:
        return
    if pulse_infos.ndim == 1:
        pulse_infos = pulse_infos.reshape(1,pulse_infos.size)

    pn = 0
    for pulse_info in pulse_infos:
        pn += 1

        pulse_time = pulse_info[2]
        pulse_sigma = pulse_info[1]
        pulse_width = pulse_info[4]
        pulse_sample = int(pulse_info[3])
        if pulse_width > 0:
            print 'Begin producing pulse number ' + str(pn) +' spectrogram figure.'


            pulse_fig_title = fig_title + '\n time = ' + str(pulse_time) + ' s ' + \
                ' , index = ' + str(int(pulse_info[3])) + \
                ' , sigma = ' + str(pulse_sigma) + \
                ' , width = ' + str(pulse_width)
            fig, axarr = setup_fig(labels, pulse_fig_title)


            chtags = ['ch1', 'ch2', 'ch3', 'ch4']
            for t in range(len(chtags)):
                ch = chtags[t]

                tag = "S1"
                fn = filepath + "_"  + ch + "_" + tag + "_new.fil"
                ext_fn = fn[:-4]+'_'+str(pulse_sample)+".ext"
                fftfreq, freq, fftdata = spectrogram(ext_fn)

                #fftdata = fftdata[:,154:156]
                #freq = freq[154:156]

#                if ch == 'ch4':
#                    print freq
#                    quit()

                fftfs = fftdata.sum(axis=1)
                fftfs = fftfs - fftfs.min()
                fftfs = np.ma.masked_less_equal(fftfs, 0)
                fftfs = 10*np.log10(fftfs)
                axarr[t][0].plot(fftfreq, fftfs)

                fftdata = fftdata - fftdata.min()
                fftdata = np.ma.masked_less_equal(fftdata, 0)
                fftdata = 10*np.log10(fftdata)


                flat_data = fftdata.flatten()
                vmin    = np.min(flat_data)
                vmax    = np.max(flat_data)
                vmean   = np.mean(flat_data)
                vmedian = np.median(flat_data)
                vstd    = np.std(flat_data)
                zmin = vmean - 3*vstd
                zmax = vmean + 3*vstd

                im = axarr[t][1].imshow(fftdata.T, aspect = 'auto',
                                   extent=[min(fftfreq),max(fftfreq), min(freq)/1000.0, max(freq)/1000.0],
                                   vmin=zmin, vmax=zmax)


                rect = axarr[t][1].get_position()
                rect.x0 = rect.x1 + 0.005
                rect.x1 = rect.x0 + 0.015
                cax = plt.axes([rect.x0,rect.y0,rect.x1-rect.x0,rect.y1-rect.y0])
                cbar = fig.colorbar(im, cax=cax, ticks = [zmin, zmax], format=sfmt, extend='both')#, aspect=5, extend='both', pad=0.0, use_gridspec=False, anchor=(-10.0,0.0))
                cbar.ax.tick_params(labelsize=10)
                cbar.ax.yaxis.get_offset_text().set_fontsize(10)



                tag = "SK"
                fn = filepath + "_"  + ch + "_" + tag + "_new.fil"
                ext_fn = fn[:-4]+'_'+str(pulse_sample)+".ext"
                fftfreq, freq, fftdata = spectrogram(ext_fn)

                fftfs = fftdata.sum(axis=1)
                fftfs = fftfs - fftfs.min()
                fftfs = np.ma.masked_less_equal(fftfs, 0)
                fftfs = 10*np.log10(fftfs)
                axarr[t][2].plot(fftfreq, fftfs)

                fftdata = fftdata - fftdata.min()
                fftdata = np.ma.masked_less_equal(fftdata, 0)
                fftdata = 10*np.log10(fftdata)

                flat_data = fftdata.flatten()
                vmin    = np.min(flat_data)
                vmax    = np.max(flat_data)
                vmean   = np.mean(flat_data)
                vmedian = np.median(flat_data)
                vstd    = np.std(flat_data)
                zmin = vmean - 3*vstd
                zmax = vmean + 3*vstd

                im = axarr[t][3].imshow(fftdata.T, aspect = 'auto',
                                   extent=[min(fftfreq),max(fftfreq), min(freq)/1000.0, max(freq)/1000.0],
                                   #norm=colors.SymLogNorm(linthresh=1.5,vmin=zmin, vmax=zmax))
                                   vmin=zmin, vmax=zmax)

                rect = axarr[t][3].get_position()
                rect.x0 = rect.x1 + 0.005
                rect.x1 = rect.x0 + 0.015
                cax = plt.axes([rect.x0,rect.y0,rect.x1-rect.x0,rect.y1-rect.y0])
                cbar = fig.colorbar(im, cax=cax, ticks = [zmin, zmax], format=sfmt, extend='both')#, aspect=5, extend='both', pad=0.0, use_gridspec=False, anchor=(-10.0,0.0))
                cbar.ax.tick_params(labelsize=10)
                cbar.ax.yaxis.get_offset_text().set_fontsize(10)
                cbar.update_ticks()



                axarr[t][3].set_xlabel('Frequency (Hz)')

            fig.savefig(fig_fn+'_'+str(pn)+'_spect.png')
            plt.close(fig)
            print 'End producing pulse number ' + str(pn) +' spectrogram figure.'

            #break
            #raw_input("Press the <ENTER> key to continue...")
        #return

#########################################
