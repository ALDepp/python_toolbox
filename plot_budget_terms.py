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

#get all data

#SST bias term
SSTbiasin = shelve.open(termspath + 'sst_daily_bias_ece-era.dat', 'r')
tos_bias_ab = SSTbiasin['tos_bias_ab']
SSTbiasin.close()

tos_bias_ab_monthly = np.asarray((np.nanmean(tos_bias_ab[:31]),np.nanmean(tos_bias_ab[31:62]),np.nanmean(tos_bias_ab[62:93]),np.nanmean(tos_bias_ab[93:])))


#Q terms
Qin = shelve.open(termspath + 'cumsum_deltaK_duetoQ.dat', 'r')
dT_ece_ab = Qin['dT_ece_ab'] #this is the temp change due to Q
dT_era_ab = Qin['dT_era_ab']
Qin.close()

dT_ece_ab_monthly = np.asarray((np.nanmean(dT_ece_ab[:31]),np.nanmean(dT_ece_ab[31:62]),np.nanmean(dT_ece_ab[62:93]),np.nanmean(dT_ece_ab[93:])))

dT_era_ab_monthly = np.asarray((np.nanmean(dT_era_ab[:31]),np.nanmean(dT_era_ab[31:62]),np.nanmean(dT_era_ab[62:93]),np.nanmean(dT_era_ab[93:])))

bias_q_ab = dT_ece_ab - dT_era_ab

bias_q_ab_monthly = np.asarray((np.nanmean(bias_q_ab[:31]),np.nanmean(bias_q_ab[31:62]),np.nanmean(bias_q_ab[62:93]),np.nanmean(bias_q_ab[93:])))

#SST dev term
SSTdevelopmentin = shelve.open(termspath + 'devSST_ab.dat', 'r')
dcum_ece_ave_ab = SSTdevelopmentin['dcum_ece_ave']
dcum_era_ave_ab = SSTdevelopmentin['dcum_era_ave']
SSTdevelopmentin.close()

dcum_ece_ave_ab_monthly = np.asarray((np.nanmean(dcum_ece_ave_ab[:31]),np.nanmean(dcum_ece_ave_ab[31:62]),np.nanmean(dcum_ece_ave_ab[62:93]),np.nanmean(dcum_ece_ave_ab[93:])))

dcum_era_ave_ab_monthly = np.asarray((np.nanmean(dcum_era_ave_ab[:31]),np.nanmean(dcum_era_ave_ab[31:62]),np.nanmean(dcum_era_ave_ab[62:93]),np.nanmean(dcum_era_ave_ab[93:])))

bias_sst_ab =  dcum_ece_ave_ab - dcum_era_ave_ab 

#horizontal advection terms
adv_ab = shelve.open('advection_terms_uv_ab.dat','r')
u_ab_ece = np.asarray(adv_ab['u_adv_terms_ab_ece'])
v_ab_ece = np.asarray(adv_ab['v_adv_terms_ab_ece'])
u_ab_ora = np.asarray(adv_ab['u_adv_terms_ab_ora'])
v_ab_ora = np.asarray(adv_ab['v_adv_terms_ab_ora'])
adv_ab.close()

uv_ece_ab = np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))
uv_ora_ab = np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))
bias_uv_ece_ora_ab = uv_ece_ab - uv_ora_ab

#upwelling term
upw_ece = shelve.open('W_below_mxl_divH_ece_ab_mesh.dat')
W_ab_ece = upw_ece['W_ab']
upw_ece.close()
contr_w_ece_ab = np.cumsum(-W_ab_ece)

upw_ora = shelve.open('W_below_mxl_divH_ora_ab_ORAmesh.dat')
W_ab_ora = upw_ora['W_ab']
upw_ora.close()
contr_w_ora_ab = np.cumsum(-W_ab_ora)

#biases and monthly means

rest_term_ece = dcum_ece_ave_ab_monthly - uv_ece_ab - (-W_ab_ece) - dT_ece_ab_monthly

rest_term_era = dcum_era_ave_ab_monthly - uv_ora_ab - (-W_ab_ora) - dT_era_ab_monthly

bias_uv_q_ab = bias_uv_ece_ora_ab + bias_q_ab_monthly


bias_upwelling_ab = contr_w_ece_ab - contr_w_ora_ab

#add the biases
bias_tot_all_terms_ab = bias_upwelling_ab + bias_q_ab_monthly + bias_uv_ece_ora_ab

#get the ave monthly sst bias
bias_sst_ab_monthly = np.asarray((np.nanmean(bias_sst_ab[:31]),np.nanmean(bias_sst_ab[31:62]),np.nanmean(bias_sst_ab[62:93]),np.nanmean(bias_sst_ab[93:])))

rest_term_bias = rest_term_ece - rest_term_era

#plots to get overview
#you miss dT in the boxes for era -> you need this to validate

#plot for presentation
# add information one by one

#plt.xlim(-1,124)
#plt.plot(range(122),dcum_ece_ave_ab,label=r'EC-E',color='blue')
#plt.plot(range(122),dcum_era_ave_ab,label=r'ERA',color='red')
#plt.plot(range(123),dT_ece_ab,label=r'Q ECE',color='blue',ls='-.',lw=1.5)
#plt.plot(range(123),dT_era_ab,label=r'Q ERA',color='red',ls='-.',lw=1.5)
#plt.plot(range(123),tos_bias_ab,label='SST ECE-ERA',color='green',lw=2)
#plt.plot(range(123),(dT_ece_ab-dT_era_ab),label=r'Q ECE-ERA',color='green',lw=2,ls='-.')
#plt.legend(bbox_to_anchor=(0.3, 0.32, 0, 0.05),prop={'size':12},fancybox=True)
#plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=14)
#plt.xlabel('Days from 05-01', fontsize=14)
#plt.ylabel(r"$\Delta$ SST [K]", fontsize=14)
#for i in np.arange(-6,2.5,0.5):
  #plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
#plt.savefig('dSST_ece_ora_Q_ab_wbelowmxl_muldzdivh.pdf',dpi=500, bbox_inches='tight')
#plt.close()

#legend(first = moves to the right, second: moves to the top)
#one by one ab box
#nb1
# only dcum sst development
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,label=r'EC-E',color='blue',lw=2)
plt.plot(range(122),dcum_era_ave_ab,label=r'ERA',color='red',lw=2)
plt.plot(range(123),dT_ece_ab,label=r'Q ECE',color='blue',ls='-.',lw=1.5,alpha=0)
plt.plot(range(123),dT_era_ab,label=r'Q ERA',color='red',ls='-.',lw=1.5,alpha=0)
plt.legend(bbox_to_anchor=(0.20, 0.18, 0, 0.05),prop={'size':14},fancybox=True)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_1_dz_test.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb3
# dcum sst development and bijdrage q
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,label=r'EC-E',color='blue',lw=2)
plt.plot(range(122),dcum_era_ave_ab,label=r'ERA',color='red',lw=2)
plt.plot(range(123),dT_ece_ab,label=r'Q ECE',color='blue',ls='-.',lw=2.5)
plt.plot(range(123),dT_era_ab,label=r'Q ERA',color='red',ls='-.',lw=2.5)
plt.legend(bbox_to_anchor=(0.20, 0.18, 0, 0.05),prop={'size':14},fancybox=True)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_2_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb4 add advection information

#these are the pure advection terms, but for their contribution to dT/dt they must be multiplied with (-1)

#dcum sst and q in the background, U and V apart
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave_ab,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece_ab,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era_ab,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,76,106),np.cumsum((-u_ab_ece)),label='ECE',color='blue',marker=r'$\mathrm{U}$',s=150)
plt.scatter((13,44,76,106),np.cumsum((-v_ab_ece)),label='ECE',color='blue',marker=r'$\mathrm{V}$',s=150)
plt.scatter((16,47,76,106),np.cumsum((-u_ab_ora)),label='ORA',color='red',marker=r'$\mathrm{U}$',s=150)
plt.scatter((16,47,76,106),np.cumsum((-v_ab_ora)),label='ORA',color='red',marker=r'$\mathrm{V}$',s=150)
plt.plot((15,46,76,107),np.cumsum((-u_ab_ece)),color='blue',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-v_ab_ece)),color='blue',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-u_ab_ora)),color='red',lw=1.5,ls=':',alpha=0.4)
plt.plot((15,46,76,107),np.cumsum((-v_ab_ora)),color='red',lw=1.5,ls=':',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.18, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_3_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#nb5 with actual bijdrage u and v
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave_ab,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece_ab,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era_ab,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),uv_ece_ab,label='ECE',color='blue',marker=r'$\mathrm{UV}$',s=330)
plt.scatter((16,47,77,107),uv_ora_ab,label='ORA',color='red',marker=r'$\mathrm{UV}$',s=330)
plt.plot((13,44,74,104),(np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))),color='blue',lw=0.5,ls=':',alpha=0.4)
plt.plot((16,47,77,107),(np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))),color='red',lw=0.5,ls=':',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.22, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_4_dz.pdf',dpi=500, bbox_inches='tight')
plt.close()

#--------- put upwelling in 

plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave_ab,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece_ab,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era_ab,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),(np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))),color='blue',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.scatter((16,47,77,107),(np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))),color='red',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),(np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))),color='blue',lw=0.5,ls=':',alpha=0.2)
plt.plot((16,47,77,107),(np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))),color='red',lw=0.5,ls=':',alpha=0.2)
plt.scatter((13,44,74,104),contr_w_ece_ab,label='ECE',color='blue',marker=r'$\mathrm{W}$',s=250)
plt.plot((13,44,74,104),contr_w_ece_ab,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),contr_w_ora_ab,label='ORA',color='red',marker=r'$\mathrm{W}$',s=250)
plt.plot((13,44,74,104),contr_w_ora_ab,color='red',lw=0.5,ls='--',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.2, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,3.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_5_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

#--------- put rest term in 

#calculate rest term

plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-6,2)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(122),dcum_ece_ave_ab,color='blue',alpha=0.2,lw=2)
plt.plot(range(122),dcum_era_ave_ab,color='red',alpha=0.2,lw=2)
plt.plot(range(123),dT_ece_ab,color='blue',ls='-.',lw=2.5,alpha=0.2)
plt.plot(range(123),dT_era_ab,color='red',ls='-.',lw=2.5,alpha=0.2)
plt.scatter((13,44,74,104),(np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))),color='blue',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.scatter((16,47,77,107),(np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))),color='red',marker=r'$\mathrm{UV}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),(np.cumsum((-u_ab_ece))+np.cumsum((-v_ab_ece))),color='blue',lw=0.5,ls=':',alpha=0.2)
plt.plot((16,47,77,107),(np.cumsum((-u_ab_ora))+np.cumsum((-v_ab_ora))),color='red',lw=0.5,ls=':',alpha=0.2)
plt.scatter((13,44,74,104),contr_w_ece_ab,color='blue',marker=r'$\mathrm{W}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),contr_w_ece_ab,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),contr_w_ora_ab,color='red',marker=r'$\mathrm{W}$',s=250,alpha=0.2)
plt.plot((13,44,74,104),contr_w_ora_ab,color='red',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),rest_term_ece,label='ECE',color='blue',marker=r'$\mathrm{R^*}$',s=250)
plt.plot((13,44,74,104),rest_term_ece,color='blue',lw=0.5,ls='--',alpha=0.4)
plt.scatter((13,44,74,104),rest_term_era,label='ORA',color='red',marker=r'$\mathrm{R^*}$',s=250)
plt.plot((13,44,74,104),rest_term_era,color='red',lw=0.5,ls='--',alpha=0.4)
plt.legend(bbox_to_anchor=(0.18, 0.2, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')
  
plt.savefig('SSTbiasdev_ab_6_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

# - ---- put bias upwelling in

plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-4.5,2.5)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(123),tos_bias_ab,color='green',lw=2.5,label='SST')

plt.scatter((15,46,76,107),bias_q_ab_monthly,color='green',marker=r'$\mathrm{Q}$',s=250,alpha=0.4)
plt.plot(range(123),bias_q_ab,color='green',lw=2.5,ls='-.',label='Q',alpha=0.4)

plt.scatter((15,46,76,107),bias_uv_ece_ora_ab,label='UV',color='green',marker=r'$\mathrm{UV}$',s=250,alpha=0.4)
plt.plot((15,46,76,107),bias_uv_ece_ora_ab,color='green',lw=0.5,ls=':',alpha=0.4)

plt.scatter((15,46,76,107),bias_upwelling_ab,label='W',color='green',marker=r'$\mathrm{W}$',s=200,alpha=0.4)
plt.plot((15,46,76,107),bias_upwelling_ab,color='green',lw=0.5,ls='--',alpha=0.4)

#plt.scatter((15,46,76,107),bias_sum_ab,label='UV + Q',marker='*',s=150,color='green')

plt.scatter((15,46,76,107),bias_tot_all_terms_ab,label='UV+W+Q',marker='*',s=250,color='black',alpha=0.8)
#plt.plot((15,46,76,107),rest_term_bias,lw=0.5,ls='--',alpha=0.4,color='green')

plt.legend(bbox_to_anchor=(0.25, 0.25, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ biases in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,2.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')

plt.savefig('SSTbiasdev_ab_7_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()

#this is with the rest term, but it just doesn't look good.
plt.figure(figsize=(10,8))
plt.xlim(-1,124)
plt.ylim(-4.5,6.5)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.plot(range(123),tos_bias_ab,color='green',lw=2.5,label='SST')

plt.scatter((15,46,76,107),bias_q_ab_monthly,color='green',marker=r'$\mathrm{Q}$',s=250)
plt.plot(range(123),bias_q_ab,color='green',lw=2.5,ls='-.',label='Q')

plt.scatter((15,46,76,107),bias_uv_ece_ora_ab,label='UV',color='green',marker=r'$\mathrm{UV}$',s=250)
plt.plot((15,46,76,107),bias_uv_ece_ora_ab,color='green',lw=0.5,ls=':',alpha=0.4)

plt.scatter((15,46,76,107),bias_upwelling_ab,label='W',color='green',marker=r'$\mathrm{W}$',s=200)
plt.plot((15,46,76,107),bias_upwelling_ab,color='green',lw=0.5,ls='--',alpha=0.4)

#plt.scatter((15,46,76,107),bias_sum_ab,label='UV + Q',marker='*',s=150,color='green')

plt.scatter((15,46,76,107),rest_term_bias,label=r'$R^*$',marker=r'$\mathrm{R^*}$',s=250,color='green')
plt.plot((15,46,76,107),rest_term_bias,lw=0.5,ls='--',alpha=0.4,color='green')

plt.legend(bbox_to_anchor=(0.18, 0.24, 0, 0.05),prop={'size':14},fancybox=True,scatterpoints = 1)
plt.title(r'$\int_{}^{}{\frac{\mathrm{d}SST}{\mathrm{d}t}}$ biases in AB box',fontsize=16)
plt.xlabel('Days from 05-01', fontsize=16)
plt.ylabel(r"$\Delta$ SST [K]", fontsize=16)
for i in np.arange(-6,6.5,0.5):
  plt.axhline(i,color='gray',alpha=0.1,ls='--')

plt.savefig('SSTbiasdev_ab_8_dz_MESH.pdf',dpi=500, bbox_inches='tight')
plt.close()
