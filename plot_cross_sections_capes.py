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
#plot oras4
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

#plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)

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

# cape 2 17.3S
#ora
lat_C2 = np.where(latTA == -17.25) 
fig, axes = plt.subplots(4,3, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA,z,thetao_TA[i,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA,z,thetao_TA[i,:,lat_C2[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_1725S_150m_2000-2009_ora.png',dpi=300)
plt.close()

#ora 4 months
lat_C2 = np.where(latTA == -17.25) 
fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA,z,thetao_TA[i+4,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA,z,thetao_TA[i+4,:,lat_C2[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i+4], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_1725S_150m_2000-2009_ora_4m.png',dpi=300)
plt.close()


#ece
lat_C2_e = np.where(latTA_e == -17.25) 
fig, axes = plt.subplots(2,2, sharex=True, sharey=True,figsize=(10,15))
axes[0,0].set_ylim(0,150)
axes[0,0].set_xlim(-50,15)
axes[0,0].invert_yaxis()
for i,ax in enumerate(axes.flatten()):
  ax.contourf(lonTA_e,z_e,votemper_TA[i,:,lat_C2[0][0],:],levels=np.arange(10,30,0.25))
  cont = ax.contour(lonTA_e,z_e,votemper_TA[i,:,lat_C2[0][0],:],levels=range(19,26))
  ax.clabel(cont, inline=1, fontsize=10,fmt='%1i',)
  ax.set_title(months[i+4], fontsize=12)
 
plt.tight_layout()
plt.savefig('Thetao_1725S_150m_2000-2009_ece.png',dpi=300)
plt.close()

#now I would like to plot the difference, that means interpolating the data on the same z grid! probably on ora, because that way I'm reducing and not adding data

# from bigger to finer grid -> put ece on era grid
cross_sec1_ece = votemper_TA[:,:,lat_C1[0][0],:]

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
plt.gca().invert_yaxis()
plt.title('6S 8E PIRATA buy from May to August')

plt.plot(thetao_TA[6,:,np.where(latTA==-6)[0][0],np.where(lonTA==8.25)[0][0]],z,label='ORAS4')
plt.plot(votemper_TA[2,:,np.where(latTA==-6)[0][0],np.where(lonTA==8.25)[0][0]],z_e,label='EC-E')
plt.plot(jul_b,z_b,label='PIRATA')
plt.ylim(0,150)
plt.xlim(10,27)
plt.gca().invert_yaxis()
plt.legend(bbox_to_anchor=(1, 0.25))
plt.savefig('T_profiles_6S_8E.ps')


#cross sections along the capes
# buy 6S / 8E
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