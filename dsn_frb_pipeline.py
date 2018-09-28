import numpy as np
import datetime as dt
from sigpyproc.Readers import FilReader
import subprocess
import multiprocessing as mp
import GPUtil
import os

def unpack(fname):
        os.system("mkdir /data3/cbochenek/%s" %fname)
        cmd = "python /home/cbochenek/raw2filterbank_files/mdscc_raw2filterbank.py -c /home/cbochenek/raw2filterbank_files/mdscc_rawparam.py -i /data3/majid/mars_data/%s/ -o /data3/cbochenek/%s/" %(fname,fname)
        out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        donelist = open("/home/cbochenek/unpacked.txt",'a')
        donelist.write("%s\n" %fname)
        donelist.close()
        print out.communicate()
        return fname


def intersect(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return np.array(merged)

def extract_70m(directory, pols=[3,4]):
	print directory, pols
        scan_table = np.genfromtxt("scan.table.all",skip_header=1,dtype = int)
        mdscc_table = np.genfromtxt("mdscc.mars.date",skip_header=1,dtype=str)
        year = int(directory.split('m')[0])+2000
        day = int(directory.split('m')[1])
        date = dt.datetime(year,1,1) + dt.timedelta(day-1)
        datestr = "%s" %str(year)[-2:]
        if date.month >9:
                datestr += "%s" %(date.month)
        else:
                datestr += "0%s" %date.month
        if date.day > 9:
                datestr+="%s" %date.day
        else:
                datestr+= "0%s" %date.day
        datestr2='%s/%s/%s' %(date.month,date.day,date.year)
        print datestr
        scan_table_lines = np.where(int(datestr) == scan_table[:,3])[0]
        print datestr2
        mdscc_lines = np.where(np.logical_and(datestr2 == mdscc_table[:,-1], mdscc_table[:,6]=='DSS-63'))[0]
        if len(mdscc_lines)>1:
                start_ends = mdscc_table[mdscc_lines][:,3:5].astype(int)
                intervals = intersect(start_ends).astype(str)
        elif len(mdscc_lines) == 0:
                return None
        else:
                intervals = mdscc_table[mdscc_lines].astype(str)[:,3:5]
        print mdscc_lines
        print scan_table_lines
        print intervals
        for i in range(len(intervals)):
                ints = intervals[i]
                start = dt.datetime(year = date.year, month = date.month, day = date.day, hour = int(ints[0][:-2]),minute = int(ints[0][-2:]))
                end = dt.datetime(year = date.year, month = date.month, day = date.day, hour = int(ints[1][:-2]),minute = int(ints[1][-2:]))
                if len(scan_table_lines) > 1:
                        for j in range(len(scan_table_lines)):
                                scan_start_cand = dt.datetime(year = date.year, month = date.month, day = date.day, hour = scan_table[scan_table_lines[j]][4], minute = scan_table[scan_table_lines[j]][5], second = scan_table[scan_table_lines[j]][6])
                                if scan_start_cand == start:
                                        scan_start = scan_start_cand
                else:
                        print scan_table[scan_table_lines]
                        scan_start = dt.datetime(year = date.year, month = date.month, day = date.day, hour = scan_table[scan_table_lines][:,4], minute = scan_table[scan_table_lines][:,5], second = scan_table[scan_table_lines][:,6])
                start_dt  = (start - scan_start)
                read_time = (end-start)
                for pol in pols:
                        files = FilReader("/data3/cbochenek/%s/mdscc_ch%s_S1.fil" %(directory,pol))
                        tsamp = files.header['tsamp']
                        start_samp = (start_dt.seconds/tsamp)+1
                        read_samps = read_time.seconds/tsamp
                        out = subprocess.Popen("extract /data3/cbochenek/%s/mdscc_ch%s_S1.fil %s %s > /data3/cbochenek/%s/mdscc_ch%s_S1_70m_scan%s.fil" %(directory,pol,start_samp,read_samps,directory,pol,i), shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                        print out.communicate()
        return directory

def run_filter60hz(fname,pols = [3,4]):
        for pol in pols:
                cmd = "ls /data3/cbochenek/%s/mdscc_ch%s_S1_70m_scan* | wc -l" %(fname, pol)
                out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                loops = int(out.communicate()[0])
                for i in range(loops):
                        print "Filtering %s poln %s loop %s out of %s" %(fname, pol,i, loops)
			#if pol == 3:
                        cmd = "python /home/cbochenek/m_fb_60hzfilter.py --inputFilename /data3/cbochenek/%s/mdscc_ch%s_S1_70m_scan%s.fil --outputFilename /data3/cbochenek/%s/mdscc_ch%s_S1_scan%s_60hzfilter.fil --maxHarmonicFrequency 180.0 --outputDir /data3/cbochenek/%s/ --clean True" %(fname,pol,i,fname,pol,i,fname)
			print cmd
			#if pol == 4:
                        #        cmd = "python /home/cbochenek/m_fb_60hzfilter.py --inputFilename /data/cbochenek/%s/mdscc_ch%s_S1_70m_scan%s.fil --outputFilename /data3/cbochenek/%s/mdscc_ch%s_S1_scan%s_60hzfilter.fil --maxHarmonicFrequency 180.0 --outputDir /data3/cbochenek/%s/ --clean True" %(fname,pol,i,fname,pol,i,fname)
                        out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                        print out.communicate()
                        print "Done filtering %s poln %s loop %s out of %s" %(fname, pol,i,loops)
                donelist = open("/home/cbochenek/filtered.txt",'a')
                donelist.write("%s %s\n" %(fname,pol))
                donelist.close()
        return (fname, pols)

def run_zap(fname, chans = [3,4]):
        cmd = "ls /data3/cbochenek/%s/mdscc_ch%s_S1_scan* | wc -l" %(fname, chans[0])
        out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        loops = int(out.communicate()[0])
        print loops
        print "Filtering RFI for %s on channels %s for %s scans" %(fname,chans,loops)
        for i in range(loops):
                cmd = "python /home/cbochenek/zap_RFI.py -f %s -s %s " %(fname,i)
                for chan in chans:
                        cmd += "-p %s " %chan
                out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                msg = out.communicate()
                print msg
        donelist = open("/home/cbochenek/cleaned.txt",'a')
        donelist.write("%s\n" %fname)
        donelist.close()
        print "Done filtering RFI for %s on channels %s" %(fname,chans)
        return fname

def get_gpu():
        cmd = "nvidia-smi"
        out = subprocess.Popen(cmd,shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        msg = out.communicate()
        ID = msg[0].split('\n')[-3].find(" 0 ")
        if ID != -1:
                gpu = 1
        else:
                gpu = 0
        return gpu

def run_heimdall(fname):
        cmd = "ls /data3/cbochenek/%s/mdscc_70m_S1_scan*" %fname
        out = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        msg = out.communicate()
        search_files = msg[0].split('\n')[:-1]
        loops = len(search_files)
        gpu = GPUtil.getFirstAvailable()[0]
        print "Running heimdall on %s for %s loops on gpu %s" %(fname, loops, gpu)
        i = 0
        for f in search_files:
                print "loop %s" %i
                out = subprocess.Popen(["/home/cbochenek/heimdall/Applications/heimdall", "-gpu_id", "%s" %gpu, "-dm", "0", "10000", "-dm_tol", "1.005", "-f", "%s" %(f)],stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                print out.communicate()
                i+=1
        print "Done running heimdall on %s for %s loops on gpu %s" %(fname, loops, gpu)
        return fname

def run_pipeline(fname, pols = [3,4]):
	#unpacked = unpack(fname)
	#extracted = extract_70m(fname)
	filtered = run_filter60hz(fname, pols)
	zapped = run_zap(fname, pols)
	searched = run_heimdall(fname)
	return fname

pool = mp.Pool(processes=16)
tracks = ["17m100","17m101","17m103","17m105","17m107","17m118","17m122","17m146","17m147"]
#tracks = ["17m152","17m154","17m157","17m167","17m175","17m181","17m182","17m183","17m195","17m196","17m345"]
#tracks = ["17m195"]#,"17m195"]#["17m181","17m182"]#,"17m183","17m195","17m196"]
results = [pool.map(run_pipeline, [tracks[i] for i in range(len(tracks))],)]
for result in results:
	print result
