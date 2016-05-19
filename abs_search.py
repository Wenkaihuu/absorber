'''
This code is designed to find the absorber candidates which have a triangle-like line profile
'''

import numpy as np
from numpy.lib.utils import safe_eval
import scipy.signal as signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pickle


parts  = ['part1','part2','part3','part4','part5','part6','part7','part8','part9']
path = '/home/p/pen/wkhu/scratch/map_result_gbt_absorb/'

sec = ''
#sec = '_secA'
#sec = '_secB'

for part in parts:
      statistic = {}
###############################   read spectra   ###########################
      #r = np.load(path+"/maps/9maps_1hr/sr_"+part+"_9parts_1hr.npy")
      #f = open(path+"/maps/9maps_1hr/sr_"+part+"_9parts_1hr.npy.meta")
      r = np.load(path+"/maps/9maps_1hr/bins/sr_bin4_"+part+"_9parts_1hr.npy")
      f = open(path+"/maps/9maps_1hr/bins/sr_"+part+"_9parts_1hr.npy.meta")

      #r = np.load(path+"data/9maps_ss/gbt_1hr_bin1"+sec+"_pointcorr_absorb_"+part+"of9/fir_1hr_bin1"+sec+"_pointcorr_absorb_"+part+"of9_clean_map_I_800.npy")
      #f = open(path+"data/9maps_ss/gbt_1hr_bin1"+sec+"_pointcorr_absorb_"+part+"of9/fir_1hr_bin1"+sec+"_pointcorr_absorb_"+part+"of9_clean_map_I_800.npy")
      
      infostring = f.readline()
      info = safe_eval(infostring)
      
      #freq_centre = info['freq_centre']/1000000000
      #freq_delta = info['freq_delta']/1000000000
      ra_centre = info['ra_centre']
      ra_delta = info['ra_delta']
      dec_centre = info['dec_centre']
      dec_delta = info['dec_delta']
      
      ra=[]
      dec=[]
      
      k=0
      while k<len(r[0][:]):
          ra.append(ra_centre+(k-len(r[0][:])//2)*ra_delta)
          k=k+1
      
      k=0
      while k<len(r[0][0][:]):
          dec.append(dec_centre+(k-len(r[0][0][:])//2)*dec_delta)
          k=k+1
      #############################   medfit   ################################
      n1=0
      while n1<len(r[0][:]):
          n2=0
          while n2<len(r[0][0][:]):
              freq_ori = np.arange(900018310.547,700030517.578,-48828.125/4)
              freq_ori /=1.e6
              temp_ori = [r[i][n1][n2] for i in range(len(r))]
      ####### throw out the >850Mhz freq which have large noise  ########        
              for i in freq_ori:
                  if np.abs(i - 850.0) < 0.048828125/4:
                      n = freq_ori.tolist().index(i)
              
              #y=signal.medfilt(y,5)
              #y=signal.medfilt(y,155)
              #medfit twice with different parameters to get the better fitted spec which has small oscillation
              fitted_y=signal.medfilt(temp_ori,15)
              minus_all = np.array(temp_ori[n:-1]) - np.array(fitted_y[n:-1])
              sigma = np.sqrt(np.var(minus_all))
              print sigma
              fitted_y=signal.medfilt(fitted_y,135)
              fitted_y = fitted_y[n:-1]
              x = np.array(freq_ori[n:-1])
              y = np.array(temp_ori[n:-1])
      
      ############################    sigma compute   ########################
      #select out one segment of the spectra which has low noise, and use the variance of this part of spectra as the variance of whole spectra to compute the sigma of the candidates

              l=len(y)-200
              while l>400 :
                  minus=[]
                  k=l
                  while k >l-401 :
                      minus.append(y[k]-fitted_y[k])
                      k=k-1
                  
                  sigma_minus=math.sqrt(np.var(minus))
                  minus = np.array(minus)
                  bol  = abs(minus)/sigma>2.5
                  if bol.sum()>0:
                      l=l-10
                  else:                         
                      break
      
      
      ##############################  find absorber and save the candidates in a dictionary  ######################
              print sigma_minus
              print part,n1,n2
              diff = (y - fitted_y).tolist()
              for i in diff[2:-2]:
                  diff_11 = diff[diff.index(i)-2]
                  diff_22 = diff[diff.index(i)+2]
                  diff_1 = diff[diff.index(i)-1]
                  diff_2 = diff[diff.index(i)+1]
                  if i<0 and np.abs(i)/sigma_minus>3.0:
                      if (diff_1<0 and np.abs(diff_1)<np.abs(i)) and (diff_2<0 and np.abs(diff_2)<np.abs(i)):
                          if (diff_11<0 and np.abs(diff_11)<np.abs(diff_1)) and (diff_22<0 and np.abs(diff_22)<np.abs(diff_2)):
                              if statistic.has_key(str(n1)+' '+str(n2)):

                                  statistic[str(n1)+' '+str(n2)].append(x[diff.index(i)])
                              else:
                                  statistic[str(n1)+' '+str(n2)]=[]
                                  statistic[str(n1)+' '+str(n2)].append(x[diff.index(i)])
                              #if statistic.has_key(str(x[diff.index(i)])):
                              #    statistic[str(x[diff.index(i)])].append(str(n1)+' '+str(n2))
                              #else:
                              #    statistic[str(x[diff.index(i)])] = []
                              #    statistic[str(x[diff.index(i)])].append(str(n1)+' '+str(n2))
              print n1,n2
              n2 = n2 +1
          n1 = n1 +1
      output = open('/home/p/pen/wkhu/scratch/map_result_gbt_absorb/statistic_absorber/1hr/plot_sr/candidate_bin4_'+part+'.pkl', 'wb')
      pickle.dump(statistic, output)
      output.close()
      
