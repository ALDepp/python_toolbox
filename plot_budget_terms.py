#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
from numpy import inf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import maskoceans
import os
import pickle
import shelve
import copy
import math
import glob
import matplotlib.lines as mlines

termspath='/net/bhw435/nobackup_1/users/deppenme/sens_exp/a01d/buget_terms/'
region = 'ab' # ab or itcz
#get all data
#SST bias term
SSTbiasin = shelve.open(termspath + 'sst_daily_bias_ece-era.dat', 'r')
tos_bias = SSTbiasin['tos_bias_' + region]
SSTbiasin.close()
tos_bias_monthly = np.asarray((np.nanmean(tos_bias[:31]),np.nanmean(tos_bias[31:62]),np.nanmean(tos_bias[62:93]),np.nanmean(tos_bias[93:])))

#Q terms
Qin = shelve.open(termspath + 'cumsum_deltaK_duetoQ.dat', 'r')
dT_ece = Qin['dT_ece_' + region] #this is the temp change due to Q
dT_era = Qin['dT_era_' + region]
Qin.close()

dT_ece_monthly = np.asarray((np.nanmean(dT_ece[:31]),np.nanmean(dT_ece[31:62]),np.nanmean(dT_ece[62:93]),np.nanmean(dT_ece[93:])))
dT_era_monthly = np.asarray((np.nanmean(dT_era[:31]),np.nanmean(dT_era[31:62]),np.nanmean(dT_era[62:93]),np.nanmean(dT_era[93:])))

bias_q = dT_ece - dT_era
bias_q_monthly = np.asarray((np.nanmean(bias_q[:31]),np.nanmean(bias_q[31:62]),np.nanmean(bias_q[62:93]),np.nanmean(bias_q[93:])))

#SST dev term
SSTdevelopmentin = shelve.open(termspath + 'devSST' + region + '.dat', 'r')
dcum_ece_ave = SSTdevelopmentin['dcum_ece_ave']
dcum_era_ave = SSTdevelopmentin['dcum_era_ave']
SSTdevelopmentin.close()

dcum_ece_ave_monthly = np.asarray((np.nanmean(dcum_ece_ave[:31]),np.nanmean(dcum_ece_ave[31:62]),np.nanmean(dcum_ece_ave[62:93]),np.nanmean(dcum_ece_ave[93:])))
dcum_era_ave_monthly = np.asarray((np.nanmean(dcum_era_ave[:31]),np.nanmean(dcum_era_ave[31:62]),np.nanmean(dcum_era_ave[62:93]),np.nanmean(dcum_era_ave[93:])))

bias_sst =  dcum_ece_ave - dcum_era_ave 

#horizontal advection terms
adv = shelve.open(termspath + 'advection_terms_uv' + region + '.dat','r')
u_ece = np.asarray(adv['u_adv_terms_' + region + '_ece'])
v_ece = np.asarray(adv['v_adv_terms_' + region + '_ece'])
u_ora = np.asarray(adv['u_adv_terms_' + region + '_ora'])
v_ora = np.asarray(adv['v_adv_terms_' + region + '_ora'])
adv.close()

uv_ece = np.cumsum((-u_ece))+np.cumsum((-v_ece))
uv_ora = np.cumsum((-u_ra))+np.cumsum((-v_ora))
bias_uv_ece_ora = uv_ece - uv_ora

#upwelling term
upw_ece = shelve.open(termspath + 'W_below_mxl_divH_ece_' + region + '_mesh.dat')
W_ece = upw_ece['W_' + region]
upw_ece.close()
contr_w_ece = np.cumsum(-W_ece)

upw_ora = shelve.open(termspath + 'W_below_mxl_divH_ora_' + region + '_ORAmesh.dat')
W_ora = upw_ora['W_' + region]
upw_ora.close()
contr_w_ora = np.cumsum(-W_ora)

#biases and monthly means
rest_term_ece = dcum_ece_ave_monthly - uv_ece - (-W_ece) - dT_ece_monthly
rest_term_era = dcum_era_ave_monthly - uv_ora - (-W_ora) - dT_era_monthly
bias_uv_q = bias_uv_ece_ora + bias_q_monthly
bias_upwelling = contr_w_ece - contr_w_ora

#add the biases
bias_tot_all_terms = bias_upwelling + bias_q_monthly + bias_uv_ece_ora

#get the ave monthly sst bias
bias_sst_monthly = np.asarray((np.nanmean(bias_sst[:31]),np.nanmean(bias_sst[31:62]),np.nanmean(bias_sst[62:93]),np.nanmean(bias_sst[93:])))
rest_term_bias = rest_term_ece - rest_term_era

#plot for presentation
#legend(first = moves to the right, second: moves to the top)
#one by one
#nb1
# only dcum sst development
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,label=r'EC-E',color='blue',lw=2)
plt.plot(range(122),dcum_era_ave,label=r'ERA',color='red',lw=2)
plt.plot(range(123),dT_ece,label=r'Q ECE',color='blue',ls='-.',lw=1.5,alpha=0)
plt.plot(range(123),dT_era,label=r'Q ERA',color='red',ls='-.',lw=1.5,alpha=0)
plt.legend(bbox_to_anchor=(0.20, 0.18, 0, 0.05),prop={'size':14},fancybox=True)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region.upper() + '_1_dz_test.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb2
# dcum sst development and bijdrage q
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,label=r'EC-E',color='blue',lw=2)
plt.plot(range(122),dcum_era_ave,label=r'ERA',color='red',lw=2)
plt.plot(range(123),dT_ece,label=r'Q ECE',color='blue',ls='-.',lw=2.5)
plt.plot(range(123),dT_era,label=r'Q ERA',color='red',ls='-.',lw=2.5)
plt.legend(bbox_to_anchor=(0.20, 0.18, 0, 0.05),prop={'size':14},fancybox=True)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region.upper() + '_2_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb3 add advection information
#these are the pure advection terms, but for their contribution to dT/dt they must be multiplied with (-1)
#dcum sst and q in the background, U and V apart
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_eceb,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,76,106),np.cumsum((-u_ece)),label='ECE',color='blue',marker=r'$\mathrm{U}$',s=150)
plt.scatter((13,44,76,106),np.cumsum((-v_ece)),label='ECE',color='blue',marker=r'$\mathrm{V}$',s=150)
plt.scatter((16,47,76,106),np.cumsum((-u_ora)),label='ORA',color='red',marker=r'$\mathrm{U}$',s=150)
plt.scatter((16,47,76,106),np.cumsum((-v_ora)),label='ORA',color='red',marker=r'$\mathrm{V}$',s=150)
plt.plot((15,46,76,107),np.cumsum((-u_ece)),color='blue',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-v_ece)),color='blue',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-u_ora)),color='red',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-v_ora)),color='red',lw=1.5,ls=':',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.18, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region.upper + '_3_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb4 with actual bijdrage u and v
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),uv_ece,label='ECE',color='blue',marker=r'$\mathrm{UV}$',s=330)
plt.scatter((16,47,77,107),uv_ora,label='ORA',color='red',marker=r'$\mathrm{UV}$',s=330)
plt.plot((13,44,74,104),(np.cumsum((-u_ece))+np.cumsum((-v_ece))),color='blue',lw=0.5,ls=':',alpha=0.4)
plt.plot((16,47,77,107),(np.cumsum((-u_ora))+np.cumsum((-v_ora))),color='red',lw=0.5,ls=':',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.22, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region + '_4_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb 5
#--------- put upwelling in 
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),(np.cumsum((-u_ece))+np.cumsum((-v_ece))),color='blue',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.scatter((16,47,77,107),(np.cumsum((-u_ora))+np.cumsum((-v_ora))),color='red',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),(np.cumsum((-u_ece))+np.cumsum((-v_ece))),color='blue',lw=0.5,ls=':',alpha=0.2)
plt.plot((16,47,77,107),(np.cumsum((-u_ora))+np.cumsum((-v_ora))),color='red',lw=0.5,ls=':',alpha=0.2)
plt.scatter((13,44,74,104),contr_w_ece,label='ECE',color='blue',marker=r'$\mathrm{W}$',s=250)
plt.plot((13,44,74,104),contr_w_ece,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),contr_w_ora,label='ORA',color='red',marker=r'$\mathrm{W}$',s=250)
plt.plot((13,44,74,104),contr_w_ora,color='red',lw=0.5,ls='--',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.2, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,3.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region + '_5_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

# nb 6
#--------- put rest term in 
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),(np.cumsum((-u_ece))+np.cumsum((-v_ece))),color='blue',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.scatter((16,47,77,107),(np.cumsum((-u_ora))+np.cumsum((-v_ora))),color='red',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),(np.cumsum((-u_ece))+np.cumsum((-v_ece))),color='blue',lw=0.5,ls=':',alpha=0.2)
plt.plot((16,47,77,107),(np.cumsum((-u_ora))+np.cumsum((-v_ora))),color='red',lw=0.5,ls=':',alpha=0.2)
plt.scatter((13,44,74,104),contr_w_ece,color='blue',marker=r'$\mathrm{W}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),contr_w_ece,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),contr_w_ora,color='red',marker=r'$\mathrm{W}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),contr_w_ora,color='red',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),rest_term_ece,label='ECE',color='blue',marker=r'$\mathrm{R^*}$',s=250)
plt.plot((13,44,74,104),rest_term_ece,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),rest_term_era,label='ORA',color='red',marker=r'$\mathrm{R^*}$',s=250)
plt.plot((13,44,74,104),rest_term_era,color='red',lw=0.5,ls='--',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.2, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_' + region + '_6_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

# nb 7
# - ---- put bias upwelling in
# this is with the star, not the good one
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-4.5,2.5)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(123),tos_bias,color='green',lw=2.5,label='SST')
plt.scatter((15,46,76,107),bias_q_monthly,color='green',marker=r'$\mathrm{Q}$',s=250,alpha=0.4)
plt.plot(range(123),bias_q,color='green',lw=2.5,ls='-.',label='Q',alpha=0.4)
plt.scatter((15,46,76,107),bias_uv_ece_ora,label='UV',color='green',marker=r'$\mathrm{UV}$',s=250,alpha=0.4)
plt.plot((15,46,76,107),bias_uv_ece_ora,color='green',lw=0.5,ls=':',alpha=0.4)
plt.scatter((15,46,76,107),bias_upwelling,label='W',color='green',marker=r'$\mathrm{W}$',s=200,alpha=0.4)
plt.plot((15,46,76,107),bias_upwelling,color='green',lw=0.5,ls='--',alpha=0.4)
plt.scatter((15,46,76,107),bias_tot_all_terms,label='UV+W+Q',marker='*',s=250,color='black',alpha=0.8)
plt.legend(bbox_to_anchor=(0.25, 0.25, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ biases in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')

plt.savefig('SSTbiasdev_' + region + '_7_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

# nb 8
#this is with the rest term
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-4.5,6.5)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(123),tos_bias,color='green',lw=2.5,label='SST')
plt.scatter((15,46,76,107),bias_q_monthly,color='green',marker=r'$\mathrm{Q}$',s=250)
plt.plot(range(123),bias_q,color='green',lw=2.5,ls='-.',label='Q')
plt.scatter((15,46,76,107),bias_uv_ece_ora,label='UV',color='green',marker=r'$\mathrm{UV}$',s=250)
plt.plot((15,46,76,107),bias_uv_ece_ora,color='green',lw=0.5,ls=':',alpha=0.4)
plt.scatter((15,46,76,107),bias_upwelling,label='W',color='green',marker=r'$\mathrm{W}$',s=200)
plt.plot((15,46,76,107),bias_upwelling,color='green',lw=0.5,ls='--',alpha=0.4)
plt.scatter((15,46,76,107),rest_term_bias,label=r'$R^*$',marker=r'$\mathrm{R^*}$',s=250,color='green')
plt.plot((15,46,76,107),rest_term_bias,lw=0.5,ls='--',alpha=0.4,color='green')
plt.legend(bbox_to_anchor=(0.18, 0.24, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ biases in ' + region.upper() + ' box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,6.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')

plt.savefig('SSTbiasdev_' + region + '_8_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()
