import numpy as np
from copy import deepcopy

from mdscc_raw2filterbank_utils import *

class FourChannelHeader:
    def __init__(self, header, FCH1):
        self.ch1 = deepcopy(header)
        self.ch2 = deepcopy(header)
        self.ch3 = deepcopy(header)
        self.ch4 = deepcopy(header)

        self.ch1['fch1']   = FCH1['ch1']
        self.ch2['fch1']   = FCH1['ch2']
        self.ch3['fch1']   = FCH1['ch3']
        self.ch4['fch1']   = FCH1['ch4']
        

class dataBlock:

    def __init__(self):
        self.data = np.zeros((0,0))

    def initialize(self, NCH, ns):
        self.data = np.zeros((NCH, ns))

    def initOutput(self, outpath, bname, tag):
        self.ofn = outpath + "/" + bname + tag + ".fil"

        self.ofp = open(self.ofn, "wb")

    def writeHeader(self, header):

        write_header(self.ofp, self.ofn,
                     header['fch1'], header['source_name'], header['tstart'],
                     header)

    def writeOutput(self):
        write_spectra_file(self.ofp, self.data)

    def finishOutput(self):
        self.ofp.close()

    def flip(self):
        self.data = self.data[::-1,:]

    def append(self, another):
        self.data = np.append(self.data, another.data, axis = 1)

    def refine(self, indices):
        self.data = self.data[indices,:]

    def SKCalc(self, S1, S2, f, M):
        self.data = SKCalc(S1.data,S2.data, f, M)

    def timeseries(self,another):
        self.data[0,:] = np.sum(another.data, axis=0)[:]




class FourChannelData:
    def __init__(self):
        self.ch1 = dataBlock()
        self.ch2 = dataBlock()
        self.ch3 = dataBlock()
        self.ch4 = dataBlock()

    def initialize(self, NCH, ns):
        self.ch1.initialize(NCH, ns)
        self.ch2.initialize(NCH, ns)
        self.ch3.initialize(NCH, ns)
        self.ch4.initialize(NCH, ns)

    def read(self, another, ifp, ns):
        (self.ch1.data,self.ch2.data,self.ch3.data,self.ch4.data,
         another.ch1.data,another.ch2.data,another.ch3.data,another.ch4.data) = \
         read_spectra_file(
            f_in, ns,
            self.ch1.data,self.ch2.data,self.ch3.data,self.ch4.data,
             another.ch1.data,another.ch2.data,another.ch3.data,another.ch4.data)

    def initOutput(self, outpath, bname, tag):
        self.ch1.initOutput(outpath, bname + "_ch1_", tag)
        self.ch2.initOutput(outpath, bname + "_ch2_", tag)
        self.ch3.initOutput(outpath, bname + "_ch3_", tag)
        self.ch4.initOutput(outpath, bname + "_ch4_", tag)


    def writeHeader(self, headers):
        self.ch1.writeHeader(headers.ch1)
        self.ch2.writeHeader(headers.ch2)
        self.ch3.writeHeader(headers.ch3)
        self.ch4.writeHeader(headers.ch4)

    def writeOutput(self):
        self.ch1.writeOutput()
        self.ch2.writeOutput()
        self.ch3.writeOutput()
        self.ch4.writeOutput()

    def finishOutput(self):
        self.ch1.finishOutput()
        self.ch2.finishOutput()
        self.ch3.finishOutput()
        self.ch4.finishOutput()

    def flip(self):
        self.ch1.flip()
        self.ch2.flip()
        self.ch3.flip()
        self.ch4.flip()

    def append(self, another):
        self.ch1.append(another.ch1)
        self.ch2.append(another.ch2)
        self.ch3.append(another.ch3)
        self.ch4.append(another.ch4)



    def refine(self, indices):
        self.ch1.refine(indices['ch1'])
        self.ch2.refine(indices['ch2'])
        self.ch3.refine(indices['ch3'])
        self.ch4.refine(indices['ch4'])

    def SKCalc(self, S1, S2, f_scale, M):
        self.ch1.data = SKCalc(S1.ch1.data,S2.ch1.data, f_scale['ch1'], M)
        self.ch2.data = SKCalc(S1.ch2.data,S2.ch2.data, f_scale['ch2'], M)
        self.ch3.data = SKCalc(S1.ch3.data,S2.ch3.data, f_scale['ch3'], M)
        self.ch4.data = SKCalc(S1.ch4.data,S2.ch4.data, f_scale['ch4'], M)

    def timeseries(self, another):
        self.ch1.timeseries(another.ch1)
        self.ch2.timeseries(another.ch2)
        self.ch3.timeseries(another.ch3)
        self.ch4.timeseries(another.ch4)



def read_4chan_spectra(S1,S2,ifp,ns):
    (S1.ch1.data,S1.ch2.data,S1.ch3.data,S1.ch4.data,
     S2.ch1.data,S2.ch2.data,S2.ch3.data,S2.ch4.data) = read_spectra_file(
        ifp, ns,
        S1.ch1.data,S1.ch2.data,S1.ch3.data,S1.ch4.data,
        S2.ch1.data,S2.ch2.data,S2.ch3.data,S2.ch4.data)



def SKCalc(S1,S2,f,M):
    kur = np.divide(S2,np.multiply(S1,S1)) # S2/(S1*S1)
    scaled_kur = np.divide(kur,f)# kur/f
    scaled_kur[np.isnan(scaled_kur)] = 0
    scaled_kur[np.isinf(scaled_kur)] = 0
    SK = np.multiply(((M+1.0)/(M-1.0)), # ((M+1.0)/(M-1.0))*(M*scaled_kur - 1.0)
    np.subtract(np.multiply(M,scaled_kur),1.0))
    return SK



def SK_CalcScale(S1,S2,M):
    kur = np.divide(S2,np.power(S1,2))#S2/(S1**2.0)
    kur[np.isnan(kur)] = 0.0
    kur[np.isinf(kur)] = 0.0
    f = np.multiply(((M+1.0)/2.0),np.mean(kur, axis=1, keepdims = 1))
    return f

def FourChannelData_SK_CalcScale(S1,S2,M):
    f_scale = {}
    f_scale['ch1'] = SK_CalcScale(S1.ch1.data, S2.ch1.data, M)
    f_scale['ch2'] = SK_CalcScale(S1.ch2.data, S2.ch2.data, M)
    f_scale['ch3'] = SK_CalcScale(S1.ch3.data, S2.ch3.data, M)
    f_scale['ch4'] = SK_CalcScale(S1.ch4.data, S2.ch4.data, M)
    return f_scale