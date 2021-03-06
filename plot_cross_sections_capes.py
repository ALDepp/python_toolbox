#!/usr/bin/env python

import os, pickle, math, shelve
import sys
sys.path.append( '/usr/people/deppenme/git_repo_0/python_toolbox' )
import plot_tool as pt
import funcs_for_bias as fs
import numpy as np
from mpl_toolkits.basemap import Basemap, maskoceans
from netCDF4 import Dataset
from numpy.random import uniform, seed
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt

path='/nobackup_1/users/deppenme/sens_exp/oras4/'

#oras4 data
infile = 'thetao_oras4_1m_2000-2009_grid_a01d_ymonmean_eragrid.nc'
bfile = Dataset(path + infile,'r')
thetao = bfile.variables['thetao'][:]
z = bfile.variables['deptht'][:]
lon = bfile.variables['longitude'][:]
lat = bfile.variables['latitude'][:]
bfile.close()

lon_180 = fs.shift_lon(lon)
thetao_180 = fs.shift_data(thetao,lon)

thetao_TA = fs.sellonlat_3d(lon_180,lat,-60,20,-35,35,thetao_180)
lonTA, latTA = fs.find_lon_lat(lon_180,lat,-60,20,-35,35)

#model data
mpath='/nobackup_1/users/deppenme/sens_exp/a01d/thetao/'
infile = 'votemper_a01d_1m_20000501_20090831_ensmean_grid_ERAImonmean.nc'
cfile = Dataset(mpath + infile,'r')
votemper = cfile.variables['votemper'][:]
z_e = cfile.variables['deptht'][:]
lon_e = cfile.variables['longitude'][:]
lat_e = cfile.variables['latitude'][:]
cfile.close()

lon_180_e = fs.shift_lon(lon_e)
votemper_180 = fs.shift_data(votemper,lon_e)

votemper_TA = fs.sellonlat_3d(lon_180_e,lat_e,-60,20,-35,35,votemper_180)
lonTA_e, latTA_e = fs.find_lon_lat(lon_180_e,lat_e,-60,20,-35,35)

months = ['January','Febuary','March','April','May','June','July','August','September','October','November','December']

#plot oras4 sst maps
#for single plot per file output
for i in range(12):
  map = Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0)
  map.drawcoastlines()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lonTA,latTA))
  map.fillcontinents(color='gray',lake_color='aqua')
  data = map.contourf(x,y,thetao_TA[i,0,:,:],levels=np.arange(15,30,0.25))
  plt.savefig('SST_map_35SN_2000-2009_' + months[i] + '_ora.ps',bbox_inches='tight')
  plt.close()
  
# for a panel plot
fig, axes = plt.subplots(3, 4) # rows, columns

for i,ax in enumerate(axes.flatten()):
    map=Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0,
    ax=ax)
    map.drawcoastlines()
    par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
    map.fillcontinents(color='gray',lake_color='aqua')
    ax.set_title(months[i])
    x, y = map(*np.meshgrid(lonTA,latTA))
    data = map.contourf(x,y,thetao_TA[i,0,:,:],levels=np.arange(15,30,0.25))

plt.tight_layout()
plt.savefig('sc_sst_oras4.ps',dpi=300)
plt.close()


#plot ece
# for individual plots
# remember that there are only 4 months in the file
for i in range(4):
  map = Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0)
  map.drawcoastlines()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lonTA_e,latTA_e))
  map.fillcontinents(color='gray',lake_color='aqua')
  data = map.contourf(x,y,votemper_TA[i,0,:,:],levels=np.arange(15,30,0.25))
  plt.savefig('SST_map_35SN_2000-2009_' + months[i+4] + '_ece.ps',bbox_inches='tight')
  plt.close()

#for a panel plot
fig, axes = plt.subplots(2,2) # rows, columns

for i,ax in enumerate(axes.flatten()):
    map=Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0,
    ax=ax)
    map.drawcoastlines()
    par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
    map.fillcontinents(color='gray',lake_color='aqua')
    ax.set_title(months[i+4])
    x, y = map(*np.meshgrid(lonTA_e,latTA_e))
    data = map.contourf(x,y,votemper_TA[i,0,:,:],levels=np.arange(15,30,0.25))

plt.tight_layout()
plt.savefig('sc_sst_ece.ps',dpi=300)
plt.close()
  
#plot the bias 
#first layer
#single plots
for i in range(4):
  map = Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0)
  map.drawcoastlines()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lonTA_e,latTA_e))
  map.fillcontinents(color='gray',lake_color='aqua')
  data = map.contourf(x,y,(-thetao_TA[i+4,0,:,:] + votemper_TA[i,0,:,:]),levels=np.arange(-3,3,0.125))
  plt.savefig('SST_bias_35SN_2000-2009_' + months[i+4] + '_ece-ora.ps',bbox_inches='tight')
  plt.close()  

#first layer
#panel plot
fig, axes = plt.subplots(2,2) # rows, columns

for i,ax in enumerate(axes.flatten()):
    map=Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0,
    ax=ax)
    map.drawcoastlines()
    par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
    map.fillcontinents(color='gray',lake_color='aqua')
    ax.set_title(months[i])
    x, y = map(*np.meshgrid(lonTA,latTA))
    data = map.contourf(x,y,(-thetao_TA[i+4,0,:,:] + votemper_TA[i,0,:,:]),levels=np.arange(-3,3,0.125))

plt.tight_layout()
plt.savefig('sc_sst_bias_ece-oras4.ps',dpi=300)
plt.close()

  
#plot the bias 
#first couple of layers - 45m
for i in range(4):
  map = Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0)
  map.drawcoastlines()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lonTA_e,latTA_e))
  map.fillcontinents(color='gray',lake_color='aqua')
  data = map.contourf(x,y,(- np.nanmean(thetao_TA[i+4,:5,:,:],axis=0) + np.nanmean(votemper_TA[i,:6,:,:],axis=0)),levels=np.arange(-3,3,0.125))
  cbar = map.colorbar(data)
  plt.savefig('0-50m_bias_35SN_2000-2009_' + months[i+4] + '_ece-ora.ps',bbox_inches='tight')
  plt.close()  


#couple of layers
#panel plot
fig, axes = plt.subplots(2,2) # rows, columns

for i,ax in enumerate(axes.flatten()):
    map=Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0,
    ax=ax)
    map.drawcoastlines()
    par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
    map.fillcontinents(color='gray',lake_color='aqua')
    ax.set_title(months[i])
    x, y = map(*np.meshgrid(lonTA,latTA))
    data = map.contourf(x,y,(- np.nanmean(thetao_TA[i+4,:5,:,:],axis=0) + np.nanmean(votemper_TA[i,:6,:,:],axis=0)),levels=np.arange(-3,3,0.125))

plt.tight_layout()
plt.savefig('sc_50m_bias_ece-oras4.ps',dpi=300)
plt.close()

#only subsurface
#panel plot
fig, axes = plt.subplots(2,2) # rows, columns

for i,ax in enumerate(axes.flatten()):
    map=Basemap(projection='cea',
    resolution='i',
    llcrnrlat=-35,urcrnrlat=35,
    llcrnrlon=-60,urcrnrlon=20,lat_ts=0,
    ax=ax)
    map.drawcoastlines()
    par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    mer = map.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])
    map.fillcontinents(color='gray',lake_color='aqua')
    ax.set_title(months[i])
    x, y = map(*np.meshgrid(lonTA,latTA))
    data = map.contourf(x,y,(- np.nanmean(thetao_TA[i+4,1:6,:,:],axis=0) + np.nanmean(votemper_TA[i,2:7,:,:],axis=0)),levels=np.arange(-3,3,0.125))

plt.tight_layout()
plt.savefig('sc_15-50m_bias_ece-oras4.ps',dpi=300)
plt.close()

#cross sections along the capes
# cape 1 0.9S
#ora
lat_C1 = np.where(latTA == -0.75) 
fig, axes = plt.subplots(4,3, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA,z,thetao_TA[i,:,lat_C1[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA,z,thetao_TA[i,:,lat_C1[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_075S_150m_2000-2009_ora.png',dpi=300)
plt.close()

#ora 4 months
lat_C1 = np.where(latTA == -0.75) 
fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA,z,thetao_TA[i+4,:,lat_C1[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA,z,thetao_TA[i+4,:,lat_C1[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i+4], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_075S_150m_2000-2009_ora_4m.png',dpi=300)
plt.close()


#ece
lat_C1_e = np.where(latTA_e == -0.75) 
fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA_e,z_e,votemper_TA[i,:,lat_C1[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA_e,z_e,votemper_TA[i,:,lat_C1[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i+4], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_075S_150m_2000-2009_ece.png',dpi=300)
plt.close()

## cape 2 17.3S ! CAPE TWO IS NOT SO INTERESTING
##ora
#lat_C2 = np.where(latTA == -17.25) 
#fig, axes = plt.subplots(4,3, sharex=True, sharey=True,figsize=(10,15))
#axes[0,0].set_ylim(0,150)
#axes[0,0].set_xlim(-50,15)
#axes[0,0].invert_yaxis()
#for i,ax in enumerate(axes.flatten()):
  #ax.contourf(lonTA,z,thetao_TA[i,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  #cont = ax.contour(lonTA,z,thetao_TA[i,:,lat_C2[0][0],:],levels=range(19,26))
  #ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  #ax.set_title(months[i], fontsize=12)
 
#plt.tight_layout()
#plt.savefig('Thetao_1725S_150m_2000-2009_ora.png',dpi=300)
#plt.close()

##ora 4 months
#lat_C2 = np.where(latTA == -17.25) 
#fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
#axes[0,0].set_ylim(0,150)
#axes[0,0].set_xlim(-50,15)
#axes[0,0].invert_yaxis()
#for i,ax in enumerate(axes.flatten()):
  #ax.contourf(lonTA,z,thetao_TA[i+4,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  #cont = ax.contour(lonTA,z,thetao_TA[i+4,:,lat_C2[0][0],:],levels=range(19,26))
  #ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  #ax.set_title(months[i+4], fontsize=12)
 
#plt.tight_layout()
#plt.savefig('Thetao_1725S_150m_2000-2009_ora_4m.png',dpi=300)
#plt.close()


##ece
#lat_C2_e = np.where(latTA_e == -17.25) 
#fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
#axes[0,0].set_ylim(0,150)
#axes[0,0].set_xlim(-50,15)
#axes[0,0].invert_yaxis()
#for i,ax in enumerate(axes.flatten()):
  #ax.contourf(lonTA_e,z_e,votemper_TA[i,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  #cont = ax.contour(lonTA_e,z_e,votemper_TA[i,:,lat_C2[0][0],:],levels=range(19,26))
  #ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  #ax.set_title(months[i+4], fontsize=12)
 
#plt.tight_layout()
#plt.savefig('Thetao_1725S_150m_2000-2009_ece.png',dpi=300)
#plt.close()

# plot subsurface differences
#fist interpolate ece on era grid (fine to coarse)
lat_C1_e = np.where(latTA_e == -0.75) 
lat_C1 = np.where(latTA == -0.75) 
cross_sec1_ece = votemper_TA[:,:,lat_C1_e[0][0],:]
cross_sec1_ora = thetao_TA[:,:,lat_C1[0][0],:]
x_ece,y_ece = np.meshgrid(lonTA_e,z_e) #(ece grid)
x,y = np.meshgrid(lonTA,z) #ora grid

cross_sec1_ece_interp = np.zeros((
		  len(cross_sec1_ece[:,0,0]),
		  len(cross_sec1_ora[0,:,0]),
		  len(cross_sec1_ora[0,0,:])))

for i in range(len(cross_sec1_ece[:,0,0])):
  cross_sec1_ece_interp[i,:,:] = griddata((x_mo.ravel(),y_mo.ravel()),
					  cross_sec1_ece[i,:,:].ravel(),
					  (x, y), method='cubic')
  
# calc diff ora and ece
diff_sec1 = np.zeros_like(cross_sec1_ece_interp)
for i in range(len(diff_sec1[:,0,0])):
  diff_sec1[i,:,:] = - cross_sec1_ora[i+4,:,:] + cross_sec1_ece_interp[i,:,:]

#plot the difference
fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  data = ax.contourf(lonTA,z,diff_sec1[i,:,:],levels=np.arange(-3,3.25,0.125))
  ax.set_title(months[i+4], fontsize=12)

plt.tight_layout()
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(data, cax=cbar_ax)
cbar.set_ticks(np.arange(-3,4,0.5))
#cbar.set_label('K')
plt.savefig('thetao_diff_075S_150m_2000-2009_ece-ora.png',dpi=300)
plt.close()

### plot where the buoys are
bpath='/nobackup_1/users/deppenme/PIRATA/buoy_8E_6S/'
bfile = 't6s8e_mon_May-Aug_2015.cdf'
bdata = Dataset(bpath + bfile, 'r')
t = bdata.variables['T_20'][:]
qual = bdata.variables['QT_5020'][:]
lat_b = bdata.variables['lat'][:]
lon_b = bdata.variables['lon'][:]
time_b = bdata.variables['time'][:]
z_b = bdata.variables['depth'][:]
bdata.close()

may_b = np.nanmean(t[:2,:,0,0],axis=0)
jun_b = np.nanmean(t[1:3,:,0,0],axis=0)
jul_b = np.nanmean(t[2:4,:,0,0],axis=0)
aug_b = np.nanmean(t[3:5,:,0,0],axis=0)


#plot the buoy data dev with time
levels = np.arange(10,30,0.25)
plt.contourf(time_b,z_b,t[:,:,0,0].transpose(),levels)
plt.ylim(0,150)
plt.ylabel('Depth')
plt.gca().invert_yaxis()
plt.title('6S 8E PIRATA buy from May to August')
plt.xlabel('Days from 01-05')
plt.savefig('buoy_6S_8E_May-Aug.ps', dpi=300)
plt.close()

#plot the buoy data dev with time -> actually 
#do the same thing for the model and oras4 so we can compare
bloc_ece = votemper_TA[:,:,np.where(latTA_e == -6)[0][0],np.where(lonTA_e==8.25)[0][0]]
levels = np.arange(10,30,0.25)
plt.contourf(range(4),z_e,bloc_ece.transpose(),levels)
plt.ylim(0,150)
plt.ylabel('Depth')
plt.gca().invert_yaxis()
plt.title('6S 8.25E ECE from May to August')
plt.xlabel('Months from 01-05')
plt.savefig('ece_6S_8E_May-Aug.ps', dpi=300)
plt.close()

bloc_ora = thetao_TA[:,:,np.where(latTA == -6)[0][0],np.where(lonTA==8.25)[0][0]]
levels = np.arange(10,30,0.25)
plt.contourf(range(4),z,bloc_ora[4:8,:].transpose(),levels)
plt.ylim(0,150)
plt.ylabel('Depth')
plt.gca().invert_yaxis()
plt.title('6S 8.25E ORAS4 from May to August')
plt.xlabel('Months from 01-05')
plt.savefig('ora_6S_8E_May-Aug.ps', dpi=300)
plt.close()

#plot the temperature profiles with depth
plt.plot(thetao_TA[6,:,np.where(latTA==-6)[0][0],np.where(lonTA==8.25)[0][0]],z,label='ORAS4')
plt.plot(votemper_TA[2,:,np.where(latTA==-6)[0][0],np.where(lonTA==8.25)[0][0]],z_e,label='EC-E')
plt.plot(jul_b,z_b,label='PIRATA')
plt.ylim(0,150)
plt.xlim(10,27)
plt.gca().invert_yaxis()
plt.legend(bbox_to_anchor=(1, 0.25))
plt.savefig('T_profiles_6S_8E.ps')

# cross section at buoy site
# buoy 6S / 8E
#ora
lat_BY = np.where(latTA == -6) 
for i in range(12):
  plt.ylim(0,150)
  plt.xlim(-50,15)
  plt.contourf(lonTA,z,thetao_TA[i,:,lat_BY[0][0],:],levels=np.arange(10,30,0.25))
  plt.gca().invert_yaxis()
  plt.title(months[i])
  plt.savefig('Thetao_6S_150m_2000-2009_' + months[i] + '_ora.ps')
  plt.close()

#ece
lat_BY_e = np.where(latTA_e == -6) 
for i in range(4):
  plt.ylim(0,150)
  plt.xlim(-50,15)
  plt.contourf(lonTA_e,z_e,votemper_TA[i,:,lat_BY_e[0][0],:],levels=np.arange(10,30,0.25))
  plt.gca().invert_yaxis()
  plt.title(months[i+4])
  plt.savefig('Thetao_6S_150m_2000-2009_' + months[i+4] + '_ece.ps')
  plt.close()
  
# plot hovmoeller diagrams
# at C1 latitude
short_months = []
for i in range(len(months)):
  short_months.append(months[i][:3])

X, Y = np.meshgrid(lonTA, range(1,13))
fig=plt.figure()
ax1 = fig.add_subplot(111)
lvls=np.arange(22,29.5,0.125)
ax1.set_xlim(-45,10)
CF = ax1.contourf(X,Y,cross_sec1_ora[:,0,:],
         levels = lvls
        )
cbar = plt.colorbar(CF, ticks=np.arange(22,31,1), format='%1i')
#cbar.set_label('K',rotation='horizontal')
plt.yticks(range(1,13),short_months)
plt.tight_layout()
plt.savefig('hovmoeller_ora.ps',dpi=300)
plt.close()

# ece data
X, Y = np.meshgrid(lonTA, range(5,9))
fig=plt.figure()
ax1 = fig.add_subplot(111)
lvls=np.arange(22,29.5,0.125)
ax1.set_xlim(-45,10)
CF = ax1.contourf(X,Y,cross_sec1_ece[:,0,:],
         levels = lvls
        )
cbar = plt.colorbar(CF, ticks=np.arange(22,31,1), format='%1i')
#cbar.set_label('K',rotation='horizontal')
plt.yticks(range(1,13),short_months)
plt.tight_layout()
plt.savefig('hovmoeller_ece.ps',dpi=300)
plt.close()

#combine
X, Y = np.meshgrid(lonTA, range(1,13))
X_e, Y_e = np.meshgrid(lonTA, range(5,9))
fig=plt.figure()
ax1 = fig.add_subplot(111)
lvls=np.arange(22,29.5,0.125)
ax1.set_xlim(-45,9)
CF = ax1.contourf(X,Y,cross_sec1_ora[:,0,:],
         levels = lvls
        )
cbar = plt.colorbar(CF, ticks=np.arange(22,31,1), format='%1i')
CF1 = ax1.contour(X,Y,cross_sec1_ora[:,0,:],
         levels = np.arange(22,30,1)
        )
plt.clabel(CF1, ticks=np.arange(22,31,1),inline=True, fmt='%.1f')
CF2 = ax1.contour(X_e,Y_e,cross_sec1_ece[:,0,:],
         levels = np.arange(22,28.5,0.5), colors='k'
        )
plt.clabel(CF2, ticks=np.arange(22,31,1),inline=True, fmt='%.1f')
plt.yticks(range(1,13),short_months)
#plt.xticks(np.arange(-50,15,10),['50W','40W','30W','20W','10W','0E','10E'])
#plt.xlabel(r'$^o$E')
#plt.tight_layout()
plt.savefig('hovmoeller_ora_ece.ps',dpi=300, bbox_inches="tight")
plt.close()
