import numpy as np
import sys
import struct
import os
import glob
import time
import datetime as dt
import matplotlib.pyplot as plt
import os.path
import string
from scipy import stats
import pandas as pd
import myutils as u

from astropy.convolution import convolve, Box1DKernel

class Filterbank:
    'common base class for all filterbanks'
    '''
    source_name = ""
    telescope_id = 0
    machine_id = 0
    data_type = 1
    fch1 = 0.0
    foff = 0.0
    nchans = 0
    obits = 0
    mjd = 0.0
    tsamp = 0.0
    nifs = 1
    barycentric = 0
    rawdatafile = ""
    az_start = 0
    za_start = 0
    src_raj = 0.0
    src_dej = 0.0
    tstart = 0.0
    period = 0.0
    ibeam = 0
    nbeams = 1
    nbits = 32
    nbins = 0
    nifs = 1
    refdm = 0.0
    '''
    
    def __init__(self, fn,readAllFlag=1):

        self.rawdatafile = fn
        totalbytes = 4
        f = open(fn, "rb")
        self.ifp = f
        string, nbytes = get_string(f)
 
        tot=0

        if (string != "HEADER_START"):
            print "error: bad type"
            exit(0)

        totalbytes=totalbytes + nbytes

        self.debug = 0

        self.source_name = ""
        self.rawdatafile = ""
        self.nifs = 0
        self.nbits = 0
        self.tstart = 0.0
        self.tsamp = 0.0
        self.fch1 = 0.0
        self.foff = 0.0
        self.nchans = 0
        self.telescope_id = 0
        self.machine_id = 0
        self.data_type = 0
        self.nbeams = 0
        self.barycentric = 0

        self.totalbytes = 0

        expecting_rawdatafile=0
        # default src raj,dej Vela Pulsar
        self.setRaDec("08:35:20.6", "-45:10:34.8")

        while (1):
            string, nbytes = get_string(f)
            #print "nbytes string = ", nbytes, string
            if (string=="HEADER_END"):
                break
            if (string=="rawdatafile"):
                expecting_rawdatafile=1
            elif (string=="source_name"):
                expecting_source_name=1
            elif (string=="FREQUENCY_START"):
                expecting_frequency_table=1
                channel_index=0
            elif (string=="FREQUENCY_END"):
                expecting_frequency_table=0
            elif (string=="az_start"):
                self.az_start, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="za_start"):
                self.za_start, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="src_raj"):
                self.src_raj, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="src_dej"):
                self.src_dej, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="tstart"):
                self.tstart, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="tsamp"):
                self.tsamp, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="period"): 
                period, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="fch1"):
                self.fch1, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="fchannel"): 
                frequency_table, = struct.unpack('f', f.read(8))
                totalbytes=totalbytes + 8
                fch1=0.0
                foff=0.0 # set to 0.0 to signify that a table is in use
            elif (string=="foff"):
                self.foff, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (string=="nchans"):
                self.nchans, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="telescope_id"):
                self.telescope_id, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="machine_id"):
                self.machine_id, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="data_type"):
                self.data_type, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="ibeam"):
                self.ibeam, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="nbeams"):
                self.nbeams, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="nbits"):
                self.nbits, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="barycentric"):
                barycentric, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="pulsarcentric"):
                pulsarcentric, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="nbins"):
                nbins, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="nsamples"): 
                itmp, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="nifs"): 
                self.nifs, = struct.unpack('i', f.read(4))
                totalbytes=totalbytes + 4
            elif (string=="npuls"): 
                npuls, = struct.unpack('d', f.read(8)) # this is actually a long int
                totalbytes=totalbytes + 8
            elif (string=="refdm"): 
                refdm, = struct.unpack('d', f.read(8))
                totalbytes=totalbytes + 8
            elif (expecting_rawdatafile): 
                self.rawdatafile = string
                expecting_rawdatafile=0
            elif (expecting_source_name):
                self.source_name = string
                expecting_source_name=0
            else: 
                print "read_header - unknown parameter: %s" %(string)
                exit(1)
            totalbytes=totalbytes + nbytes

        totalbytes=totalbytes + len("HEADER_END")
        self.totalbytes = totalbytes
        #print "totalbytes", self.totalbytes
        #print "tot", tot 

        fileSizeInBytes = os.path.getsize(fn)
        dataSizeInBytes = fileSizeInBytes - self.totalbytes

        ns = dataSizeInBytes/4
        self.nsamps = ns/self.nchans
        print "file size = %d, ns = %d" % (fileSizeInBytes, ns)
        print "nsamps = %d" % (self.nsamps)
        print "nchans = %d" % (self.nchans)

        if (readAllFlag):
            fmt = get_fmt_string(ns)
            mydata = f.read(ns*4)
            self.rawData = np.fromstring(mydata, dtype='<f4')
            print "size of mydata ", len(mydata)
            print "size of rawData ", len(self.rawData)

            self.rawData = self.rawData.reshape(self.nsamps, self.nchans)
            self.bp_median = np.median(self.rawData, axis=0)
            self.bp = np.mean(self.rawData, axis=0)
            self.ts = np.mean(self.rawData, axis=1)
            self.dur = self.nsamps * self.tsamp
            self.corrData = np.zeros((self.nsamps, self.nchans))

            # capture statistics of data points
            flat_data = self.rawData.flatten()
            self.vmin    = np.min(flat_data)
            self.vmax    = np.max(flat_data)
            self.vmean   = np.mean(flat_data)
            self.vmedian = np.median(flat_data)
            self.vstd    = np.std(flat_data)

            print "stats: min,max,mean,median,rms = ", self.vmin,self.vmax,self.vmean,self.vmedian,self.vstd

    def printHeader(self):
        myname = "print_header"
        print "rawdatafile  = ", self.rawdatafile
        print "source_name  = ", self.source_name
        print "telescope_id = ", self.telescope_id
        print "machine_id   = ", self.machine_id
        print "data_type    = ", self.data_type
        print "nifs         = ", self.nifs
        print "nbits        = ", self.nbits
        print "tstart       = ", self.tstart
        print "tsamp        = ", self.tsamp
        print "fch1         = ", self.fch1
        print "foff         = ", self.foff
        print "nchans       = ", self.nchans
        print "nbeams       = ", self.nbeams
        print "barycentric  = ", self.barycentric

    def set_bp_ts(self):
        self.bp_median = np.median(self.rawData, axis=0)
        self.bp = np.mean(self.rawData, axis=0)
        self.ts = np.mean(self.rawData, axis=1)
        self.dur = self.nsamps * self.tsamp

    def char_bp(self, win):

        if (win <= 0):
            print "error: set_bp: illegal value for win"
            exit(0)

        nbp_threshold = 3
        bad_ch = []

        # you have to have read some data already

        # first estimate the corrected bp value
        bp1 = self.bp
        bptmp = bp1[0:win]
        bp2 = np.hstack((bptmp,bp1))
        bp_smooth = pd.rolling_median(bp2,win)[win:]
        self.corrbp = bp_smooth

        # obtain a metric by which i can identify bad channels
        corrbptmp = self.corrbp[0:win]
        corrbp2 = np.hstack((corrbptmp,self.corrbp))
        bp_std = pd.rolling_std(corrbp2,win)[win:]
        self.bp_std = bp_std
        bp_delta = (self.bp - self.corrbp)
        bp_delta_norm = (bp_delta - np.mean(bp_delta))/np.std(bp_delta)

        # identify bad channels
        self.nbp_delta = bp_delta_norm
        ii = np.where(abs(self.nbp_delta) > nbp_threshold)
        ii_badch = ii[0]
        print "set_bp: nbad, bad channels: ", len(ii_badch), ii_badch
        for ch in ii_badch:
            bad_ch.append(ch)

        self.bad_ch = bad_ch

        # set good_ch
        self.good_ch = np.ones(self.nchans, dtype=np.int)
        for ch in self.bad_ch:
            self.good_ch[ch] = 0


    def test_bad_fpga_chans(self, bp):
        d  = bp - np.roll(bp, 1)
        dn = (d - np.mean(d))/np.std(d)
        ii = np.where(abs(dn) > 1)
        ii_array = ii[0]
        n = len(ii_array)
        print "bad_fpga_chans: n = ", n
        print "bad_fpga_chans: ", ii_array
        if (n>1):
            return ii_array[1]
        else:
            return 0

    # find bad channels resulting from the FPGA FFT engine
    # the characteristic of these channels is a sharp drop off in counts to almost 0
    def bad_fpga_chans(self, bp, plotit=0):
        DELTA_THRESHOLD=-2.5
        bad_ch = []

        xi = np.arange(self.nchans)
        y = bp
        slope, intercept, r_val, p_val, std_err = stats.linregress(xi, y)
        line = slope*xi + intercept
        print "bad_fpga_chans: fit result: r,p,rms: ", r_val, p_val, std_err

        d = y - line
        dn = (d - np.mean(d))/np.std(d)

        if plotit:
            plt.figure()
            plt.subplot(211)
            plt.title('bad_fpga_chans')
            plt.plot(xi, line, 'r-', xi, y,'o')
            plt.grid()
            plt.subplot(212)
            plt.plot(xi,dn)
            plt.axhline(y=DELTA_THRESHOLD, xmin=0, xmax=self.nchans - 1, c='r',linewidth=0.5, zorder=0)
            plt.grid()
            plt.show()
        
        ii = np.where(dn < DELTA_THRESHOLD)
        ii_badch = ii[0]
        for ch in ii_badch:
            bad_ch.append(ch)

        return bad_ch

    # find bad channels - outliers
    def bad_chans(self, bp, plotit=0):
        DELTA_THRESHOLD=2.0
        bad_ch = []

        xi = np.arange(self.nchans)
        y = bp
        slope, intercept, r_val, p_val, std_err = stats.linregress(xi, y)
        line = slope*xi + intercept
        print "bad_chans: fit result: r,p,rms: ", r_val, p_val, std_err

        d = y - line
        dn = (d - np.mean(d))/np.std(d)

        if plotit:
            plt.figure()
            plt.subplot(211)
            plt.title('bad_chans')
            plt.plot(xi, line, 'r-', xi, y,'o')
            plt.grid()
            plt.subplot(212)
            plt.plot(xi,dn)
            plt.axhline(y=DELTA_THRESHOLD, xmin=0, xmax=self.nchans - 1, c='r',linewidth=0.5, zorder=0)
            plt.axhline(y=-1*DELTA_THRESHOLD, xmin=0, xmax=self.nchans - 1, c='r',linewidth=0.5, zorder=0)
            vmin = np.min(dn); vmax = np.max(dn)
            if vmin > -1*DELTA_THRESHOLD:
                vmin = -1*DELTA_THRESHOLD - 1
            if vmax < DELTA_THRESHOLD:
                vmax = DELTA_THRESHOLD + 1
            plt.ylim(vmin - (0.1*np.abs(vmin)), vmax + (0.1*np.abs(vmax)))
            plt.grid()
            plt.show()
        
        ii = np.where(abs(dn) > DELTA_THRESHOLD)
        ii_badch = ii[0]
        for ch in ii_badch:
            bad_ch.append(ch)

        return bad_ch

    def bpCorr7(self, chunkData): 

        # it looks like it is better to just correct the bp with the 
        # mean values, i.e., self.bp
        # on the other hand, to identify bad channels, i need the 
        # rolling medium with a window of at least 5
        for i in np.arange(self.nchans):
            self.corrData[:,i] = chunkData[:,i] / self.bp[i]
            #self.corrData[:,i] = chunkData[:,i] / self.corrbp[i]

        # finally zero the bad channels
        for ch in self.bad_ch:
            self.corrData[:,ch] = 0

        # capture statistics of good data points
        good_data = self.corrData[:, self.good_ch]
        flat_good_data = good_data.flatten()
        self.vmin    = np.min(flat_good_data)
        self.vmax    = np.max(flat_good_data)
        self.vmean   = np.mean(flat_good_data)
        self.vmedian = np.median(flat_good_data)
        self.vstd    = np.std(flat_good_data)

        if self.debug:
            print "bpCorr7: stats: min,max,mean,median,rms = ", self.vmin,self.vmax,self.vmean,self.vmedian,self.vstd

    def bad_fpga_chans2(self, bp):
        DELTA_THRESHOLD=2.0
        
        fillVal = np.median(bp)

        nbad = len(bp)
        i = 0
        bad_ch = []

        while (nbad > 1):
            d  = bp - np.roll(bp, 1)
            dn = (d - np.mean(d))/np.std(d)
            ii = np.where(abs(dn) > DELTA_THRESHOLD)
            ii_array = ii[0]
            n = len(ii_array)
            print "bad_fpga_chans: ", ii_array
            i = i+1

            if (i > 10):
                print "bad_fpga_chans: too many bad channels; gave up"
                exit(1)

            if (n>1):
                ch = ii_array[1]  # pick off first bad channel, skipping zeroth channel
                print "bad_fpga_chans: ch fillVal ", ch, fillVal
                bp[ch] = fillVal
                bad_ch.append(ch)
            else:
                break
        return bad_ch

    def plotEvent(self,t,sm,type):
        plt.figure()
        s = int(t/self.tsamp)
        dt = 50.0e-3 # 50 ms
        ds = int(dt/self.tsamp)
        print "event: ", t,dt
        print "event: ", s-ds, s+ds
        tstart = t - dt
        tend   = t + dt

        datain  = self.rawData[s-ds:s+ds]
        dataout = datain
        for i in range(sm):
            dataout = u.avg2d2(datain)
            datain = dataout

        if type==0:
            plt.imshow(dataout.T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans])
            #plt.imshow(self.rawData[s-ds:s+ds].T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans])
            #plt.imshow(self.rawData[s-ds:s+ds].T, aspect='auto',extent=[tstart*1000,tend*1000, 0, self.nchans])
            #plt.imshow(self.rawData[s-ds:s+ds].T, aspect='auto', extent=[0,self.dur * 1000, 0, self.nchans])
        elif type==1:
            zmin = self.vmean - 3*self.vstd
            zmax = self.vmean + 3*self.vstd
            #plt.imshow(self.rawData[s-ds:s+ds].T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans],vmin=zmin,vmax=zmax)
            plt.imshow(dataout.T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans],vmin=zmin,vmax=zmax)
        elif type==2:
            zmin = self.vmean - 5*self.vstd
            zmax = self.vmean + 5*self.vstd
            plt.imshow(dataout.T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans],vmin=zmin,vmax=zmax)
        elif type==3:
            zmin = self.vmin
            zmax = self.vmax
            plt.imshow(dataout.T, aspect='auto',extent=[-dt*1000,dt*1000, 0, self.nchans],vmin=zmin,vmax=zmax)
        plt.show()

    def plotRaw(self,type):
        plt.figure()
        if type==0:
            plt.imshow(self.rawData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans])
        elif type==1:
            zmin = self.vmean - 3*self.vstd
            zmax = self.vmean + 3*self.vstd
            plt.imshow(self.rawData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        elif type==2:
            zmin = self.vmean - 5*self.vstd
            zmax = self.vmean + 5*self.vstd
            plt.imshow(self.rawData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        elif type==3:
            zmin = self.vmin
            zmax = self.vmax
            plt.imshow(self.rawData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        plt.show()

    def plotCorr(self,type):
        plt.figure()
        if type==0:
            plt.imshow(self.corrData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans])
        elif type==1:
            zmin = self.vmean - 3*self.vstd
            zmax = self.vmean + 3*self.vstd
            plt.imshow(self.corrData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        elif type==2:
            zmin = self.vmean - 5*self.vstd
            zmax = self.vmean + 5*self.vstd
            plt.imshow(self.corrData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        elif type==3:
            zmin = self.vmin
            zmax = self.vmax
            plt.imshow(self.corrData.T, aspect='auto', extent=[0,self.dur, 0, self.nchans],
                       vmin=zmin, vmax=zmax)
        plt.show()

    def writeFil(self, fn):
        writeHeader(fn)
        #writeData(self, fn)

    def setRaDec(self, ra, dec):
        els = dec.split(":")
        
        sgn = int(els[0])
        if sgn < 0:
            sgn=-1
        else:
            sgn=1
                
        decVal = np.abs(int(els[0])*10000) + int(els[1])*100 + float(els[2])
        decVal = sgn*decVal

        els = ra.split(":")
                
        raVal = np.abs(int(els[0])*10000) + int(els[1])*100 + float(els[2])

        self.src_raj = raVal
        self.src_dej = decVal

    def extractChans(self, loCh, hiCh):
        self.rawData = self.rawData[:, loCh:hiCh]
        self.nchans = hiCh - loCh
        
        self.bp = np.sum(self.rawData, axis=0)/self.nsamps
        self.ts = np.sum(self.rawData, axis=1)/self.nchans
        self.dur = self.nsamps * self.tsamp
        self.corrData = np.zeros((self.nsamps, self.nchans))


    def writeHeader(self, ofn):
        ofp = open(ofn, "wb")
        send_string(ofp, "HEADER_START")

        send_string(ofp, "rawdatafile")
        send_string(ofp, ofn)

        send_string(ofp, "source_name")
        send_string(ofp, self.source_name)

        send_int(ofp, "telescope_id",self.telescope_id)
        send_int(ofp, "machine_id",self.machine_id)
        send_int(ofp, "data_type",self.data_type) # filterbank
        send_double(ofp, "fch1",self.fch1)
        send_double(ofp, "foff",self.foff)
        send_int(ofp, "nchans",self.nchans)
        send_int(ofp, "nbits",self.nbits) # obits?
        send_double (ofp, "tstart",self.tstart) 
        send_double (ofp, "src_raj",self.src_raj) 
        send_double (ofp, "src_dej",self.src_dej) 
        send_double(ofp, "tsamp",self.tsamp)
        send_int(ofp, "nifs",self.nifs)
        send_int(ofp, "barycentric",self.barycentric)

        send_string(ofp, "HEADER_END")

        ofp.close()

         
    def writeData(self,ofn,flag):
        ofp = open(ofn, "a+b")
        if flag:
            data = self.corrData.flatten()
        else:
            data = self.rawData.flatten()

        n = len(data)
        # write smaller chucks at a time
        chunkSize = 1000000
        nchunks = n / chunkSize
        print "nchunks %d" % (nchunks)

        for i in range(nchunks):
            if not i % 100:
                print "chunk %d" % (i+1)
            istart = i*chunkSize
            iend   = istart + chunkSize
            d = data[istart:iend]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

        remainderChunk = n % chunkSize
        istart = nchunks*chunkSize
        #istart = (nchunks-1)*chunkSize
        iend   = istart + remainderChunk
        print "write: istart iend", istart, iend
        d = data[istart:iend]
        fmt = '<' + str(len(d)) + 'f'
        obuf = struct.pack(fmt, *d)
        ofp.write(obuf)

        ofp.close()

    def write2(self,ofn,data1,data2):
        ofp = open(ofn, "a+b")

        [ns, nch] = np.shape(data1)

        for i in range(ns):
            d = data1[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data2[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

        ofp.close()

    def write6(self,ofn,data1,data2,data3,data4,data5,data6):
        ofp = open(ofn, "a+b")

        [ns, nch] = np.shape(data1)

        for i in range(ns):
            d = data1[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data2[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data3[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data4[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data5[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

            d = data6[i,:]
            fmt = '<' + str(len(d)) + 'f'
            obuf = struct.pack(fmt, *d)
            ofp.write(obuf)

        ofp.close()

    def readDataChunk(self, nsamps):
        ns = nsamps * self.nchans
        if self.debug:
            print "ns = ", ns
        fmt = get_fmt_string(ns)
        mydata = self.ifp.read(ns*4)
        if self.debug:
            print "shape of mydata = ", np.shape(mydata)
        data = np.fromstring(mydata, dtype='<f4')
        if self.debug:
            print "shape of data = ", np.shape(data)
        data = data.reshape(nsamps, self.nchans)
        return data

    def genStatsChunk(self, data):
        self.bp_median = np.median(data, axis=0)
        self.bp = np.mean(data, axis=0)
        self.ts = np.mean(data, axis=1)
        self.corrData = np.zeros(np.shape(data))

        # capture statistics of data points
        flat_data = data.flatten()
        self.vmin    = np.min(flat_data)
        self.vmax    = np.max(flat_data)
        self.vmean   = np.mean(flat_data)
        self.vmedian = np.median(flat_data)
        self.vstd    = np.std(flat_data)

        print "stats: min,max,mean,median,rms = ", self.vmin,self.vmax,self.vmean,self.vmedian,self.vstd

    def writeDataChunk(self, data,ofn):
        ofp = open(ofn, "a+b")
        data = self.corrData.flatten()
        fmt = '<' + str(len(data)) + 'f'
        obuf = struct.pack(fmt, *data)
        ofp.write(obuf)
        ofp.close()

    def readChunk(self, dur):
        ns = self.nsamps
        nc = self.nchans
        dt = self.tsamp

        nsInChunk = int(np.floor(dur / self.tsamp))
        print "Will read %f second of data at a time: nsInChunk %d" % (dur, nsInChunk)
        nchunks = int(np.floor(ns / nsInChunk))
        print "nchunks = ", nchunks

        data = self.readDataChunk(nsInChunk)
        plt.imshow(data.T, aspect='auto', extent=[0,dur, 0, self.nchans])
        return data

    def procFile(self, ifn, ofn,ra,dec,dt,bp_win):
        self.setRaDec(ra, dec)
        self.tsamp = dt
        ns = self.nsamps
        nc = self.nchans
        dt = self.tsamp
        self.writeHeader(ofn)

        chunkDur = 60 # process 60s chunks at a time
        nsInChunk = int(np.floor(chunkDur / self.tsamp))
        print "Will process %f second of data at a time: nsInChunk %d" % (chunkDur, nsInChunk)
        #nsInChunk = 500000
        nchunks = int(np.floor(ns / nsInChunk))
        print "nchunks = ", nchunks
        #print "nchunks %d" % (nchunks)
        first = 1

        for i in range(nchunks):
            #if not i % 100:
            #    print "chunk %d" % (i+1)
            if self.debug:
                print "processing chunk %d ..." % (i+1)
            # read data
            self.chunkData = self.readDataChunk(nsInChunk)
            self.rawData = self.chunkData
            
            # obtain bp info and correction
            if (first):
                self.set_bp_ts()
                self.char_bp(bp_win)

            # generate statistics
            self.genStatsChunk(self.chunkData)
            # apply bandpass corr
            self.bpCorr7(self.chunkData)
            # write out corrected data
            self.writeDataChunk(self.chunkData, ofn)
            first = 0

        rChunk = ns % nsInChunk
        print "rChunk %d" % (rChunk)
        if (rChunk):
            # read data
            self.chunkData = self.readDataChunk(rChunk)
            self.rawData = self.chunkData

            # obtain bp info and correction
            if (first):
                self.set_bp_ts()
                self.char_bp(bp_win)

            # generate statistics
            self.genStatsChunk(self.chunkData)
            # apply bandpass corr
            self.bpCorr7(self.chunkData)
            # write out corrected data
            self.writeDataChunk(self.chunkData, ofn)
            first = 0

        
def get_string(f):
    n, = struct.unpack('i',f.read(4))
    # struct.unpack return tuple. how to get data type in native format
    # either use a , like above, or x = n[0]
    string = f.read(n)
    nbytes = 4 + n
    return string, nbytes

def get_fmt_string(n):
    fmt = '<' + str(n) + 'f'
    #print fmt
    return fmt

def send_string(ofp, str):
    ofp.write(struct.pack("<i",len(str)))
    ofp.write(str)

def send_int(ofp, str, value):
    send_string(ofp, str)
    ofp.write(struct.pack("<i",value))

def send_double(ofp, str, value):
    send_string(ofp, str)
    ofp.write(struct.pack("<d",value))


def corrFil(ifn, ofn, kernel, blank,ra,dec):
    f1 = Filterbank(ifn)
    f1.setRaDec(ra, dec)
    Filterbank.writeCorr(f1, ofn, kernel, blank)

#def main(ifn, ofn,ra,dec,dt, bp_win):
#    f1 = Filterbank(ifn, 0)  # create object, but don't read the data yet
#    dt = dt / 1.0e6
#    print "sample time is = %f" % (dt)
#    Filterbank.procFile(f1, ifn, ofn,ra,dec,dt, bp_win)

# usage:
# python ~/src/py/filterbank_class.py "infile" "outfile" 0 0 "05:31:47.0" "+33:05:46.5" dt bp_win >& log.sl &

if __name__ == '__main__':
    from sys import argv
    main(argv[1], argv[2],argv[3], argv[4],float(argv[5]), int(argv[6]))
    #main(argv[1], argv[2], int(argv[3]), int(argv[4]),argv[5], argv[6],float(argv[7]))
