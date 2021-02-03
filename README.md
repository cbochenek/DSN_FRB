# DSN_FRB

This repo contains an analysis pipeline that searches for FRBs from the Deep Space Network.

dsn_frb_pipeline.py - This is the script that ties everything together. It runs the steps of the pipeline on multiple observations in parallel. The steps that are run are:  
1) Unpacking the raw data into filterbank format. (extract_mars_data.py)   
2) Select the appropriate data from the 70m telescope (as opposed to data from other telescopes). (Read and interpret the information in mdscc.mars.date and scan.table.all)  
3) Search for the pestilence that is 60 Hz noise and filter it out, if present. (m_fb_60hzfilter.py)  
4) Correct for the passband, remove RFI, and baseline the data.  (zap_RFI.py)  
5) Search the data for FRBs with heimdall.  

make_plots.py - Reads the output of heimdall and makes a PDF file of all the diagnostic plots needed to decide if a candidate is real
