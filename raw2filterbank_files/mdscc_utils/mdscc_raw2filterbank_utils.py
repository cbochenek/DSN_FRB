import numpy as np
import struct
## Read scan file into list dictionary
def read_scan_file(scan_file_path):
    fs= open(scan_file_path)

    scan = []
    for line in fs:
        line = line.strip()
        columns = line.split()
        source_info = {}
        source_info['name'] = columns[0]
        source_info['duration_min'] = columns[1]
        source_info['duration_sec'] = columns[2]
        source_info['date'] = columns[3]
        source_info['start_hour'] = columns[4]
        source_info['start_min'] = columns[5]
        source_info['start_sec'] = columns[6]
        source_info['tag'] = columns[7]
        scan.append(source_info)

    return scan

## Read source scan start time from scan file
def get_tsart(source_info):

    source_date = str(source_info['date'])
#    source_date = source_info['date']
    source_start_hour = int(source_info['start_hour'])
    source_start_min = int(source_info['start_min'])
    source_start_sec = int(source_info['start_sec'])


    year = int(source_date[0:2])+2000
    month  = int(source_date[2:4])
    day  = int(source_date[4:6])
    hh   = source_start_hour
    mm   = source_start_min
    ss   = source_start_sec


    a = (14 - month)/12
    y = year + 4800 - a
    m = month + 12 * a - 3

    jdn = day + (153 * m + 2)/5 + 365 * y +y/4 - y/100 + y/400 - 32045
    mjd = jdn - 2400000.5

    mjd_frac = (hh + (mm/60.0) + (ss/3600.0))/24.0
    mjd = int(mjd) + mjd_frac

    tstart = mjd

    return tstart

# Read data packet function
def read_spectra_file(fp,nspectra, mom2_band1, mom2_band2, mom2_band3,mom2_band4,
mom4_band1, mom4_band2, mom4_band3, mom4_band4):

    Nch = 8
    Ndat = 1024
    tmp = np.fromfile(fp,dtype=np.dtype('>H'),count = (Nch*Ndat+8)*nspectra);
    tmp = tmp.reshape(2*nspectra,(Nch*Ndat/2)+4)
    tmp = tmp[:,4:]
    tmp = tmp.reshape(Nch*Ndat*nspectra)
    tmp = tmp.reshape(nspectra,2*Ndat,Nch/2)

    mom2_band1[0::2,:] = tmp[:,0:512,0].transpose()
    mom2_band1[1::2,:] = tmp[:,512:1024,0].transpose()

    mom4_band1[0::2,:] = tmp[:,0:512,1].transpose()
    mom4_band1[1::2,:] = tmp[:,512:1024,1].transpose()

    mom2_band2[0::2,:] = tmp[:,0:512,2].transpose()
    mom2_band2[1::2,:] = tmp[:,512:1024,2].transpose()

    mom4_band2[0::2,:] = tmp[:,0:512,3].transpose()
    mom4_band2[1::2,:] = tmp[:,512:1024,3].transpose()

    mom2_band3[0::2,:] = tmp[:,0+Ndat:512+Ndat,0].transpose()
    mom2_band3[1::2,:] = tmp[:,512+Ndat:1024+Ndat,0].transpose()

    mom4_band3[0::2,:] = tmp[:,0+Ndat:512+Ndat,1].transpose()
    mom4_band3[1::2,:] = tmp[:,512+Ndat:1024+Ndat,1].transpose()

    mom2_band4[0::2,:] = tmp[:,0+Ndat:512+Ndat,2].transpose()
    mom2_band4[1::2,:] = tmp[:,512+Ndat:1024+Ndat,2].transpose()

    mom4_band4[0::2,:] = tmp[:,0+Ndat:512+Ndat,3].transpose()
    mom4_band4[1::2,:] = tmp[:,512+Ndat:1024+Ndat,3].transpose()

    return  (mom2_band1, mom2_band2, mom2_band3, mom2_band4,
             mom4_band1, mom4_band2, mom4_band3, mom4_band4)


def write_spectra_file(f_out, data):
    datatmp = data.transpose().reshape(data.size)
    fmt = '<' + str(datatmp.size) + 'f'
    obuf = struct.pack(fmt, *datatmp)
    f_out.write(obuf)

## Write header to Output file
def write_header(ofp, ofn, fch1, source_name, tstart, header):

    telescope_id = header['telescope_id']
    machine_id = header['machine_id']
    data_type = header['data_type']
    foff = header['foff']
    nchans = header['nchans']
    obits = header['obits']
    tsamp = header['tsamp']
    nifs = header['nifs']
    barycentric = header['barycentric']


    write_string(ofp, "HEADER_START")

#    write_string(ofp, "rawdatafile")
#    write_string(ofp, ofn)

    write_string(ofp, "source_name")
    write_string(ofp, source_name)

    write_int(ofp, "telescope_id",telescope_id)
    write_int(ofp, "machine_id",machine_id)
    write_int(ofp, "data_type",data_type) # filterbank
    write_double(ofp, "fch1",fch1)
    write_double(ofp, "foff",foff)
    write_int(ofp, "nchans",nchans)
    write_int(ofp, "nbits",obits)
    write_double (ofp, "tstart",tstart)
    write_double(ofp, "tsamp",tsamp)
    write_int(ofp, "nifs",nifs)
    write_int(ofp, "barycentric",barycentric)

    write_string(ofp, "HEADER_END")

def write_string(ofp, str):
    ofp.write(struct.pack("<i",len(str)))
    ofp.write(str)

def write_int(ofp, str, value):
    write_string(ofp, str)
    ofp.write(struct.pack("<i",value))

def write_double(ofp, str, value):
    write_string(ofp, str)
    ofp.write(struct.pack("<d",value))

def write_spectra(f_out, data):
    fmt = '<' + str(len(data)) + 'f'
    obuf = struct.pack(fmt, *data)
    f_out.write(obuf)
