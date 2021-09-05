#!/usr/bin/python

# this code is written for python >= 3.6 and require matplotlib and numpy
# usage: python rdatafile.py inputfile outputfile
# inputfile = file name of the CADI data input file
# outputfile = file name of the resulting plot
#
# outputfile should have extension .png to produce png-files or .pdf to produce pdf-files.
# supported formats eps, pdf, pgf, png, ps, raw, rgba, svg, svgz (depends on version of matplotlib)
#
# the code is written based on IDL code provided by Chris Meek
# IDL code originally written by Ian Grant and modified by the same and JWM
#

import sys
import struct
import datetime
from time import strptime
import numpy as np
import matplotlib
import copy
# matplotlib.use('Agg')       # use the 'Agg' backend when $DISPLAY environment is not defined
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (5.12, 5.12)
import pdb;
import pandas as pd

dheight = 3.0                       # not defined in data file

import glob
Path = 'E:\Paper3-data\150318'
filenames = sorted(glob.glob('*.md4'))
          
def freading(f):
 try:
   # 1) read header information as described in the documentation p. 26-27
   site = f.read(3).decode("utf-8")
   ascii_datetime = f.read(22).decode("utf-8")
   filetype = f.read(1).decode("utf-8")
   nfreqs = struct.unpack("<H", f.read(2))[0]
   ndops = struct.unpack("<B", f.read(1))[0]
   minheight = struct.unpack("<H", f.read(2))[0]
   maxheight = struct.unpack("<H", f.read(2))[0]
   pps = struct.unpack("<B", f.read(1))[0]
   npulses_avgd = struct.unpack("<B", f.read(1))[0]
   base_thr100 = struct.unpack("<H", f.read(2))[0]
   noise_thr100 = struct.unpack("<H", f.read(2))[0]
   min_dop_forsave = struct.unpack("<B", f.read(1))[0]
   dtime = struct.unpack("<H", f.read(2))[0]
   gain_control = f.read(1).decode("utf-8")
   sig_process = f.read(1).decode("utf-8")
   noofreceivers = struct.unpack("<B", f.read(1))[0]
   spares = f.read(11).decode("utf-8")

   month = ascii_datetime[1:4]
   day = int(ascii_datetime[5:7])
   hour = int(ascii_datetime[8:10])
   minute = int(ascii_datetime[11:13])
   sec = int(ascii_datetime[14:16])
   year = int(ascii_datetime[17:21])

   month_number = strptime(month, '%b').tm_mon
   mydate = datetime.date(year, month_number, day)
   jd = mydate.toordinal() + 1721424.5
   jd0jd = datetime.date(1986, 1, 1)
   jd0 = jd0jd.toordinal() + 1721424.5

   time_header = (jd - jd0)*86400 + hour*3600 + minute*60 + sec
   time_hour = 3600 * (time_header/3600)

   max_ntimes = int(3600.0/dtime)
   max_ndopbins = 100000
   
   # 2) read all frequencies used

   freqs = [struct.unpack("<f", f.read(4))[0] for i in range(nfreqs)]

   if filetype == 'I':
      max_nfrebins = nfreqs
   else:
      max_nfrebins = min(max_ntimes * nfreqs, max_ndopbins)
      
   # 3) read rawdata

   nheights = int(maxheight / dheight + 1)

   times = []
   frebins = []
   frebins_x = []
   frebins_gain_flag = []
   frebins_noise_flag = []
   frebins_noise_power10 = []
   time_min = 0
   time_sec = 0
   timex = -1
   freqx = nfreqs - 1
   dopbinx = -1
   frebinx = -1
   iq_bytes = np.zeros((noofreceivers,2))
   dopbin_x_timex = []
   dopbin_x_freqx = []
   dopbin_x_hflag = []
   dopbin_x_dop_flag = []
   dopbin_iq = []
   hflag = 0
   var = 0
   F1fminhmin = []
   
   for num in range(6):
      
      time_min = struct.unpack("<B", f.read(1))[0]
      #while time_min != 255:
      time_sec = struct.unpack("<B", f.read(1))[0]
      flag = struct.unpack("<B", f.read(1))[0]  # gainflag
      timex += 1
      times.append(time_hour + 60 * time_min + time_sec)
      for freqx in range(nfreqs):
         noise_flag = struct.unpack("<B", f.read(1))[0] # noiseflag
         noise_power10 = struct.unpack("<H", f.read(2))[0] 
         frebinx += 1
         frebins_gain_flag.append(flag)
         frebins_noise_flag.append(noise_flag)
         frebins_noise_power10.append(noise_power10)
         flag = struct.unpack("<B", f.read(1))[0]
         while flag < 224:
            ndops_oneh = struct.unpack("<B", f.read(1))[0]
            hflag = flag
            if ndops_oneh >= 128:
               ndops_oneh = ndops_oneh - 128
               hflag = hflag + 200
            for dopx in range(ndops_oneh):
               dop_flag = struct.unpack("<B", f.read(1))[0]
               for rec in range(noofreceivers):
                  iq_bytes[rec,0] = struct.unpack("<B", f.read(1))[0]
                  iq_bytes[rec,1] = struct.unpack("<B", f.read(1))[0]
               dopbinx += 1
               dopbin_iq.append(copy.deepcopy(iq_bytes))
               dopbin_x_timex.append(timex)
               dopbin_x_freqx.append(freqx)
               dopbin_x_hflag.append(hflag)
               dopbin_x_dop_flag.append(dop_flag)
            flag = struct.unpack("<B", f.read(1))[0] # next hflag/gainflag/FF
            
      time_min = flag # next record
      #pdb.set_trace()
      
      
      frequency=[]
      for i in range(len(dopbin_x_freqx)):
          frequency.append(freqs[dopbin_x_freqx[i]]/1000000.0)

   # --- noise reduction - remove lowest and highest frequency for each height
      for index in range(320):
         vh = 30 + index
         collection = [i for i, val in enumerate(dopbin_x_hflag) if val == vh]

         if len(collection) == 1:
            del dopbin_x_hflag[collection[0]]
            del frequency[collection[0]]

         if len(collection) == 2:
            del dopbin_x_hflag[collection[0]]
            del frequency[collection[0]]
            del dopbin_x_hflag[collection[1]-1]
            del frequency[collection[1]-1]

         minimum = 13.
         maximum = 1.
         if len(collection) > 2:
            for j in range(len(collection)):
               if frequency[collection[j]] < minimum:
                  minindex = collection[j]
                  minimum = frequency[collection[j]]
               if frequency[collection[j]] > maximum:
                  maxindex = collection[j]
                  maximum = frequency[collection[j]]

            if minindex < maxindex:
               del dopbin_x_hflag[minindex]
               del frequency[minindex]
               del dopbin_x_hflag[maxindex-1]
               del frequency[maxindex-1]
            else:
               del dopbin_x_hflag[maxindex]
               del frequency[maxindex]
               del dopbin_x_hflag[minindex-1]
               del frequency[minindex-1]
      # --- end noise reduction
 
      height = list(np.array(dopbin_x_hflag) * 3)
   
      var  = 10*num
      
      title = str(month) + str(day)
      
      time = hour*3600 + var*60
      EL = []
      F1L = []
      F2L = []
      Eflag = 0
      F1flag = 0
      F2flag = 0

      for i in range(len(frequency)):
            if (height[i] > 90 and height[i] <= 140):
                EL.append((time,height[i],frequency[i]))
                Eflag = 1
            elif (height[i] > 150 and height[i] <= 400):
                F1L.append((time,height[i],frequency[i]))
                F1flag = 1
            elif (height[i] > 400 and height[i] <= 800):
                F2L.append((time,height[i],frequency[i]))
                F2flag = 1
      if (Eflag == 1):
          tempE=(min(EL, key = min))
      else: tempE = (time, 0, 0.00)
      if (F1flag == 1):
          tempF1=(min(F1L, key = min))
      else: tempF1 = (time, 0, 0.00)
      if (F2flag == 1):
          tempF2=(min(F2L, key = min))
      else: tempF2 = (time, 0, 0.00)
      
      EL = open(title + "EL" + ".txt","a+")
      EL.write("%s %d %5.2f\n" %(tempE))
      F1L = open(title + "F1L" + ".txt","a+")
      F1L.write("%s %d %5.2f\n" %(tempF1))
      F2L = open(title + "F2L" + ".txt","a+")
      F2L.write("%s %d %5.2f\n" %(tempF2))
      
      dopbinx = -1
      dopbin_x_freqx = []
      dopbin_x_hflag = []
     
   return
   EL.close()
   F1L.close()
   F2L.close()
 finally:
   f.close()

for filename in filenames:
    print(filename)
    f = open(filename,"rb")
    freading(f)
    continue


   

