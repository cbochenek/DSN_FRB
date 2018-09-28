import numpy as np
import datetime as dt
from sigpyproc.Readers import FilReader
import subprocess
import multiprocessing as mp

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

def extract_70m(directory, pols=[1,2,3,4]):
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
			files = FilReader("/data/cbochenek/%s/mdscc_ch%s_S1.fil" %(directory,pol))
			tsamp = files.header['tsamp']
			start_samp = (start_dt.seconds/tsamp)+1
			read_samps = read_time.seconds/tsamp
			out = subprocess.Popen("extract /data/cbochenek/%s/mdscc_ch%s_S1.fil %s %s > /data/cbochenek/%s/mdscc_ch%s_S1_70m_scan%s.fil" %(directory,pol,start_samp,read_samps,directory,pol,i), shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			print out.communicate()
	return directory

#directories = np.genfromtxt("unpacked.txt",dtype=str)
#print directories[2]
#extract_70m(directories[2])

pool = mp.Pool(processes=14)
#directories = np.genfromtxt("unpacked.txt",dtype=str)
directories = ["17m195"]
#results = [pool.map(extract_70m,[directories[i] for i in range(len(directories))],)]
results = extract_70m(directories[0], pols = [3,4])
outfile = open("/home/cbochenek/extracted.txt", 'a')
for result in results:
	outfile.write("%s\n" %result)
	print result
outfile.close()
