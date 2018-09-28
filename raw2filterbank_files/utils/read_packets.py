#!/usr/bin/env python

import numpy as np
import time
import matplotlib.pyplot as plt


#fp= open('./data_capture/pta_dat');
fp= open('./pta_file');

end_packet = 1;
pow_tot = np.zeros((1024,1));
pow_tot2 = np.zeros((1024,1));
pow_tot3 = np.zeros((1024,1));
pow_tot4 = np.zeros((1024,1));

kurt_tot = np.zeros((1024,1));
kurt_tot2 = np.zeros((1024,1));
kurt_tot3 = np.zeros((1024,1));
kurt_tot4 = np.zeros((1024,1));


for j in range (0, end_packet):
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
 

plt.figure()
plt.subplot(4,1,1)
plt.plot(pow_tot)
plt.subplot(4,1,2)
plt.plot(pow_tot2)
plt.subplot(4,1,3)
plt.plot(pow_tot3)
plt.subplot(4,1,4)
plt.plot(pow_tot4)

plt.figure()
plt.subplot(4,1,1)
plt.plot(kurt_tot)
plt.subplot(4,1,2)
plt.plot(kurt_tot2)
plt.subplot(4,1,3)
plt.plot(kurt_tot3)
plt.subplot(4,1,4)
plt.plot(kurt_tot4)

plt.show()
