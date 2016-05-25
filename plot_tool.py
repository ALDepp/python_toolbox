#!/usr/bin/env python

import os, pickle, math, shelve
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans

def plot_map(lon,lat,var,name,title=''):
  """
  This function plots a map of the variable you specify.
  Arguments are:
  lon, lat, variable to plot, name of file to be saved to, title (optional).
  A large range is used for the variable and the spacing.
  This applies to eg global maps. Where the range is smaller 
  plot_map_region should be used. 
  """
  mapfile = ('map_'+str(int(min(lat)))+'_'+str(int(max(lat)))+'_'
  +str(int(min(lon)))+'_'+str(int(max(lon)))+'.pkl')
  if os.path.isfile(mapfile):
    f    = open(mapfile,'r')
    map = pickle.load(f)
    f.close()
  else:
    map = Basemap(
    projection='cea',
    resolution='i',
    llcrnrlat=min(lat)+1,urcrnrlat=max(lat)-1,
    llcrnrlon=min(lon),urcrnrlon=max(lon),lat_ts=0
    )
    f = open(mapfile,'wb')
    pickle.dump(map,f)
    f.close()
  map.drawcoastlines()
  map.drawcountries()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lon,lat))
#  levels=np.arange(math.ceil(np.nanmin(var)),int(np.nanmax(var)),0.5)
  levels=np.arange(int(np.nanmax(var))-30,int(np.nanmax(var))+1,0.5)
  data = map.contourf(x,y,var,levels)
  cbar = map.colorbar(data)
  plt.title(title)
  plt.savefig(name+'.pdf',dpi=200)
  plt.close()
  return None

def plot_map_region(lon,lat,var,name,title='',lb=20,ub=30,spacing=0.25):
  """
  This function plots a map of the variable you specify.
  Arguments are:
  lon, lat, variable to plot, name of file to be saved to, title (optional).
  A large range is used for the variable and the spacing.
  This applies to eg global maps. Where the range is smaller 
  plot_map_region should be used. 
  """
  mapfile = ('map_'+str(int(min(lat)))+'_'+str(int(max(lat)))+'_'
  +str(int(min(lon)))+'_'+str(int(max(lon)))+'.pkl')
  if os.path.isfile(mapfile):
    f    = open(mapfile,'r')
    map = pickle.load(f)
    f.close()
  else:
    map = Basemap(
    projection='cea',
    resolution='i',
    llcrnrlat=min(lat)+1,urcrnrlat=max(lat)-1,
    llcrnrlon=min(lon),urcrnrlon=max(lon),lat_ts=0
    )
    f = open(mapfile,'wb')
    pickle.dump(map,f)
    f.close()
  map.drawcoastlines()
  map.drawcountries()
  par = map.drawparallels(np.arange(-80,81,5),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,2.5),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lon,lat))
  levels=np.arange(lb,ub,spacing)
  data = map.contourf(x,y,var,levels)
  cbar = map.colorbar(data)
  plt.title(title)
  plt.tight_layout()
  plt.savefig(name+'.pdf',dpi=200)
  plt.close()
  return None

def plot_bias_map(lon,lat,var1,var2,name,title='',lb=-2.5,ub=2.5,spacing=0.1):
  """
  This function plots a map of the difference between the 
  two variables you specify.
  The arguments are:
  lon, lat, var1 - var2, filename, title (optional), 
  lower and upper bound, spacing for the map (lb, ub, spacing), optional.
  """
  mapfile = ('map_'+str(int(min(lat)))+'_'+str(int(max(lat)))+'_'
  +str(int(min(lon)))+'_'+str(int(max(lon)))+'.pkl')
  if os.path.isfile(mapfile):
    f    = open(mapfile,'r')
    map = pickle.load(f)
    f.close()
  else:
    map = Basemap(
    projection='cea',
    resolution='i',
    llcrnrlat=min(lat)+1,urcrnrlat=max(lat)-1,
    llcrnrlon=min(lon),urcrnrlon=max(lon),lat_ts=0
    )
    f = open(mapfile,'wb')
    pickle.dump(map,f)
    f.close()
  map.drawcoastlines()
  map.drawcountries()
  par = map.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lon,lat))
  bias = var1 - var2
  landmask = maskoceans(x,y,np.ones(np.shape(x)),inlands=False)
  landmask.mask = ~landmask.mask
  masked_bias = bias * landmask
  levels=np.arange(lb,ub,spacing)
  data = map.contourf(x,y,bias,levels)
  cbar = map.colorbar(data)
  plt.title(title)
  plt.tight_layout()
  plt.savefig(name+'_bias.pdf',dpi=200)
  plt.close()
  return landmask

def plot_bias_map_small(lon,lat,var1,var2,name,title='',lb=-2.5,ub=2.5,spacing=0.1):
  """
  This function plots a map of the difference between the 
  two variables you specify.
  The arguments are:
  lon, lat, var1 - var2, filename, title (optional), 
  lower and upper bound, spacing for the map (lb, ub, spacing), optional.
  """
  mapfile = ('map_'+str(int(min(lat)))+'_'+str(int(max(lat)))+'_'
  +str(int(min(lon)))+'_'+str(int(max(lon)))+'.pkl')
  if os.path.isfile(mapfile):
    f    = open(mapfile,'r')
    map = pickle.load(f)
    f.close()
  else:
    map = Basemap(
    projection='cea',
    resolution='i',
    llcrnrlat=min(lat)+1,urcrnrlat=max(lat)-1,
    llcrnrlon=min(lon),urcrnrlon=max(lon),lat_ts=0
    )
    f = open(mapfile,'wb')
    pickle.dump(map,f)
    f.close()
  map.drawcoastlines()
  map.drawcountries()
  par = map.drawparallels(np.arange(-80,81,2),labels=[1,0,0,0])
  mer = map.drawmeridians(np.arange(-180,180,2),labels=[0,0,0,1])
  x, y = map(*np.meshgrid(lon,lat))
  bias = var1 - var2
  landmask = maskoceans(x,y,np.ones(np.shape(x)),inlands=False)
  landmask.mask = ~landmask.mask
  masked_bias = bias * landmask
  levels=np.arange(lb,ub,spacing)
  data = map.contourf(x,y,bias,levels)
  cbar = map.colorbar(data)
  plt.title(title)
  plt.tight_layout()
  plt.savefig(name+'_bias_s.pdf',dpi=200)
  plt.close()
  return landmask