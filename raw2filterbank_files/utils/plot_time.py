#!/usr/bin/env python

import numpy as np
import time
import matplotlib.pyplot as plt


#fp = open('./data_capture/pta_308_201610041853_001');
fp = open('./pta_file');

end_packet = 1000;
pow_tot = np.zeros((1024,1));
pow_tot2 = np.zeros((1024,1));
pow_tot3 = np.zeros((1024,1));
pow_tot4 = np.zeros((1024,1));

kurt_tot = np.zeros((1024,1));
kurt_tot2 = np.zeros((1024,1));
kurt_tot3 = np.zeros((1024,1));
kurt_tot4 = np.zeros((1024,1));

pow_sum1 = np.zeros((end_packet,1));
pow_sum2 = np.zeros((end_packet,1));
pow_sum3 = np.zeros((end_packet,1));
pow_sum4 = np.zeros((end_packet,1));

for k in range (0, end_packet):
    #seq_no = fread(fd,1,'ubit62')
    #mux_status = fread(fd,1,'ubit2')
    discard = np.fromfile(fp, dtype=np.dtype('>Q'), count = 1)

    # 4096 x 16bit numbers.

    # due to misalignment of data (same as in BRAM)
    # the 1 & 3 samples are for input 1,
    # the 2 & 4 samples are for input 2.
    j=0;
    for i in range (0, 512):      
        pow_tot[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1); 
        kurt_tot[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        pow_tot2[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        kurt_tot2[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);    
   	j=j+2;
  
    j=1;
    for i in range (0, 512):      
        pow_tot[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);         
        kurt_tot[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        pow_tot2[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        kurt_tot2[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        j=j+2;
        

    discard = np.fromfile(fp, dtype=np.dtype('<Q'), count = 1)
    j=0;
    for i in range (0, 512):      
        pow_tot3[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);         
        kurt_tot3[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        pow_tot4[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        kurt_tot4[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        j=j+2;
    
    j=1;
    for i in range (0, 512):      
        pow_tot3[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);        
        kurt_tot3[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        pow_tot4[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);
        kurt_tot4[j] = np.fromfile(fp,dtype=np.dtype('>H'),count = 1);        
        j=j+2;
 
    pow_sum1[k] = np.sum(pow_tot);
    pow_sum2[k] = np.sum(pow_tot2);
    pow_sum3[k] = np.sum(pow_tot3[100:600]);
    pow_sum4[k] = np.sum(pow_tot4[200:400]) + np.sum(pow_tot4[600:800]);
 

plt.figure()
plt.subplot(4,1,1)
plt.plot(pow_sum1);
plt.subplot(4,1,2)
plt.plot(pow_sum2);
plt.subplot(4,1,3)
plt.plot(pow_sum3);
plt.subplot(4,1,4)
plt.plot(pow_sum4);

plt.show();
