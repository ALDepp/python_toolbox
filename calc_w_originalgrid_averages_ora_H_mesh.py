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

# output dimensions: 292 lats and 362 lons
R = 6371229 #earth radius according to NEMO

def sellonlat(lon,lat,lon_min,lon_max,lat_min,lat_max,data):
  """
  Finds the lat and lon indices for the box you want.
  Assums that lat goes from 90 to -90
  and lon from -180 to 180.
  #Then selects the data for the box you want.
  """
  lomin = np.min(np.where(lon > lon_min))
  lomax = np.max(np.where(lon < lon_max))
  lamax = np.max(np.where(lat > lat_min))
  lamin = np.min(np.where(lat < lat_max))
  cut_data = data[:,lamin:lamax+1,lomin:lomax+1]
  return cut_data

flowpath = '/net/bhw435/nobackup_1/users/deppenme/sens_exp/oras4/'
meshpath = '/nobackup_1/users/deppenme/mesh/'

#get u data 
ufile = 'uo_oras4_1m_2000-2009_grid_U_ymonmean_MJJA_a01d.nc'
uin = Dataset(flowpath + ufile,'r')
nlon_ora_u = uin.variables['nav_lon'][:]
nlat_ora_u = uin.variables['nav_lat'][:]
u_ora = uin.variables['uo'][:]
uin.close()

#get v data 
vfile = 'vo_oras4_1m_2000-2009_grid_V_ymonmean_MJJA_a01d.nc'
vin = Dataset(flowpath + vfile,'r')
nlon_ora_v = vin.variables['nav_lon'][:]
nlat_ora_v = vin.variables['nav_lat'][:]
v_ora = vin.variables['vo'][:]
vin.close()

#get T data
Tfile = 'thetao_oras4_1m_2000-2009_grid_a01d_ymonmean_MJJA.nc'
Tin = Dataset(flowpath + Tfile,'r')
dept_ora_t = Tin.variables['deptht'][:]
nlon_ora_t = Tin.variables['nav_lon'][:]
nlat_ora_t = Tin.variables['nav_lat'][:]
temp_ora = Tin.variables['thetao'][:]
Tin.close()


# get the mesh
#this is calculated with the 42 layer grid - use that one

#mfile = 'mesh_zgr_42.nc' #this is the wrong grid, this is the ece grid
mfile = 'mesh_mask_ORAS4_Hao.nc'
m_in = Dataset(meshpath + mfile, 'r')
e3t = m_in.variables['e3t'][0,:,:,:] # vertical T metric partial steps ! this is the one I want
e3w = m_in.variables['e3w'][0,:,:,:]
m_in.close()

# dw/dz = - (du/dx + dv/dz)

#get dx and dy
#for dx you need lats
latt_263 = nlat_ora_t[:,263]

#def dx(lon) function depending on lat
#this function assumes a distance of 1degree always
#which is true in my region
def dx_func(phi):
  dx = 2 * R * np.pi * np.cos(np.deg2rad(phi)) / 360
  return dx

#for dy you need the interpolated diff_lat of the grid boxes
#first lats at v levels (interpolate)
latv = latt_263[:-1]/2 +  latt_263[1:]/2 
diff_lat_v = latv[1:] - latv[:-1]

#define dy function depending on diff lat
def dy_func(diff_lat):
  dy =  2 * R * np.pi / (360) *(diff_lat)
  return dy

##get d(uT) - > loop over lat
#latt_ab = latt_263[113:129]
#latt_itcz = latt_263[153:161]

#grid_box_lats_ab = grid_box_lats[113:129]
#grid_box_lats_itcz = grid_box_lats[153:161]

#diff_latt_ab = diff_latt_263[113:129]
#diff_latt_itcz = diff_latt_263[153:161]

# get du/dx
du_x = u_ora[:,:,:,1:] - u_ora[:,:,:,:-1]

#loop over lat to get diff dx's
dudx = np.zeros_like(du_x)
for lat in range(len(du_x[0,0,:,0])):
  dudx[:,:,lat,:] = du_x[:,:,lat,:] / dx_func(latt_263[lat])

# get dv/dy
dv_y = v_ora[:,:,1:,:] - v_ora[:,:,:-1,:]

#loop over lat to get diff lats - 
#the result will only be true in your region!!
dvdy = np.zeros_like(dv_y[:,:,:-1,:])
for lat in range(len(dv_y[0,0,:-1,0])):
  dvdy[:,:,lat,:] = dv_y[:,:,lat,:] / dy_func(diff_lat_v[lat])

#you have dvdy and dudx but not for T0 ! only from T1 onwards.

dwdz = - (dudx[:,:,:-2,:] + dvdy[:,:,:,:-1])
#now you have dwdz from T1 onwards. 

#multiply with dz

dw = np.zeros_like(dwdz)
for mon in range(len(dwdz[:,0,0,0])):
  for lat in range(len(dwdz[0,0,:,0])):
    for lon in range(len(dwdz[0,0,0,:])):
      dw[mon,:,lat,lon] = dwdz[mon,:,lat,lon] * e3t[:,lat,lon]

#now accumulate dw from the bottom, you will have all dws up to the top.
w = np.zeros_like(dw)
for mon in range(len(dw[:,0,0,0])):
  for lat in range(len(dw[0,0,:,0])):
    for lon in range(len(dw[0,0,0,:])):
      w[mon,:,lat,lon] = np.cumsum(dw[mon,::-1,lat,lon])[::-1]

#now I got w. I have to calculate dT/dz
dT_z = temp_ora[:,1:,:,:] - temp_ora[:,:-1,:,:]

dTdz = np.zeros_like(dT_z)
for mon in range(len(dT_z[:,0,0,0])):
  for lat in range(len(dT_z[0,0,:,0])):
    for lon in range(len(dT_z[0,0,0,:])):
      dTdz[mon,:,lat,lon] = dT_z[mon,:,lat,lon] / e3t[:-1,lat,lon]

#got get w*dT/dz (it is not possible everywhere, you don't have dT at k=0
#and you don't have w at i,j = 0
#so use  your dTdz from lon,lat = 1 and not zero

wdTdz = w[:,1:,:,:] * dTdz[:,:,1:-1,1:]

#now weigh the layer that you care about in the box that you care about
#remember all the indices have to be moved to the left.

#i want the transport below? -> into the mxl
idx_mxl = np.zeros_like(dTdz[:,0,1:-1,1:])
right_shape_dTdz = dTdz[:,:,1:-1,1:]
for mon in range(len(right_shape_dTdz[:,0,0,0])):
  for lat in range(len(right_shape_dTdz[0,0,:,0])):
    for lon in range(len(right_shape_dTdz[0,0,0,:])):
      try:
	idx_mxl[mon,lat,lon] = np.min(np.where(right_shape_dTdz[mon,:,lat,lon]<-0.1))
      except ValueError:
	idx_mxl[mon,lat,lon] = np.nan

#the depths I have now are BELOW the mixed layer, not the bottom of the mixed layer
#I think i have to use the bottom of, so when I use this index I have to say
# idx - 1 !!

#now it's time to accumulate the dwTdz multiplied with their volumes
#also sum over the volumes so that you can divide by it later. 
#maybe it's smart to first get your box, you won't encounter any nans there
#idx already uses the right grid, just move the indices 
#the real coordinates are
#lat : itcz = 152:161 ; ab =  113:129
#lon : itcz = 263:280 ; ab = 292:299
#but remember that you have to move them by -1 because your array starts at 1

idx_mxl_itcz = idx_mxl[:,151:160,262:279]
idx_mxl_ab = idx_mxl[:,112:128,291:298]

wdTdz_itcz = wdTdz[:,:,151:160,262:279]
wdTdz_ab = wdTdz[:,:,112:128,291:298]

#now i have to get the right lats and diff_lats again -> one to the right from wdTdz and idx_mxl indices
latt_itcz = latt_263[152:161]
latt_ab = latt_263[113:129]

latv_diff_itcz = diff_lat_v[152:161]
latv_diff_ab = diff_lat_v[113:129]

e3t_itcz = e3t[:,152:161,263:280]
e3t_ab = e3t[:,113:129,292:299]

def vbox_func(dz,dlat,lat):
  vbox = dz * dy_func(dlat) * dx_func(lat)
  return vbox

# first of all: use below ! the mixed layer, not at the mixed layer
# mult contr with dz and divide by h because you're warming the whole layer

w_wdTdz_boxjes_itcz = np.zeros_like(idx_mxl_itcz)
volume = np.zeros_like(w_wdTdz_boxjes_itcz)
for mon in range(len(idx_mxl_itcz[:,0,0])):
  for lat in range(len(idx_mxl_itcz[0,:,0])):
    for lon in range(len(idx_mxl_itcz[0,0,:])):
      try:
	idx_depth = (idx_mxl_itcz[mon,lat,lon])
	w_wdTdz_boxjes_itcz[mon,lat,lon] = wdTdz_itcz[mon,idx_depth,lat,lon] * vbox_func(e3t_itcz[idx_depth,lat,lon],latv_diff_itcz[lat],latt_itcz[lat])
	volume[mon,lat,lon] = vbox_func(dept_ora_t[idx_depth],latv_diff_itcz[lat],latt_itcz[lat])
      except IndexError:
	w_wdTdz_boxjes_itcz[mon,lat,lon] = np.nan
	volume[mon,lat,lon] = np.nan

wdTdz_itcz_may = np.nansum(w_wdTdz_boxjes_itcz[0,:,:],axis=(0,1)) / np.nansum(volume[0,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_itcz_jun = np.nansum(w_wdTdz_boxjes_itcz[1,:,:],axis=(0,1)) / np.nansum(volume[1,:,:],axis=(0,1)) * 24 * 3600 * 30

wdTdz_itcz_jul = np.nansum(w_wdTdz_boxjes_itcz[2,:,:],axis=(0,1)) / np.nansum(volume[2,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_itcz_aug = np.nansum(w_wdTdz_boxjes_itcz[3,:,:],axis=(0,1)) / np.nansum(volume[3,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_itcz_weighted = np.asarray((wdTdz_itcz_may,wdTdz_itcz_jun,wdTdz_itcz_jul,wdTdz_itcz_aug))

# weighted upwelling in itcz at bottom mixed layer 
#[ 0.28233164,  0.46530265,  0.34003224,  0.38273516] # i dont know where these come from - very weird
#new
#array([ 0.82667283,  0.723076  ,  0.71597365,  0.33715435])
# weighted upwelling in itcz below mixed layer
#[ 0.42826538,  1.37962839,  0.96703019,  0.46361545]
#new
#array([ 1.73702512,  1.75346754,  1.25142722,  0.96643908])

outt = shelve.open('W_below_mxl_divH_ora_itcz_ORAmesh.dat') #(idx -1)
outt['W_itcz'] = wdTdz_itcz_weighted
outt.close()

#---------------------------------------------------------

#same approach for ab region
# div by H

w_wdTdz_boxjes_ab_b = np.zeros_like(idx_mxl_ab)
volume_b = np.zeros_like(w_wdTdz_boxjes_ab_b)
for mon in range(len(idx_mxl_ab[:,0,0])):
  for lat in range(len(idx_mxl_ab[0,:,0])):
    for lon in range(len(idx_mxl_ab[0,0,:])):
      try:
	idx_depth = (idx_mxl_ab[mon,lat,lon])
	w_wdTdz_boxjes_ab_b[mon,lat,lon] = wdTdz_ab[mon,idx_depth,lat,lon] * vbox_func(e3t_ab[idx_depth,lon,lat],latv_diff_ab[lat],latt_ab[lat])
	volume_b[mon,lat,lon] = vbox_func(dept_ora_t[idx_depth],latv_diff_ab[lat],latt_ab[lat])
      except IndexError:
	w_wdTdz_boxjes_ab_b[mon,lat,lon] = np.nan
	volume_b[mon,lat,lon] = np.nan

wdTdz_ab_may_b = np.nansum(w_wdTdz_boxjes_ab_b[0,:,:],axis=(0,1)) / np.nansum(volume_b[0,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_ab_jun_b = np.nansum(w_wdTdz_boxjes_ab_b[1,:,:],axis=(0,1)) / np.nansum(volume_b[1,:,:],axis=(0,1)) * 24 * 3600 * 30

wdTdz_ab_jul_b = np.nansum(w_wdTdz_boxjes_ab_b[2,:,:],axis=(0,1)) / np.nansum(volume_b[2,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_ab_aug_b = np.nansum(w_wdTdz_boxjes_ab_b[3,:,:],axis=(0,1)) / np.nansum(volume_b[3,:,:],axis=(0,1)) * 24 * 3600 * 31

wdTdz_ab_ab_weighted_b = np.asarray((wdTdz_ab_may_b,wdTdz_ab_jun_b,wdTdz_ab_jul_b,wdTdz_ab_aug_b))

out = shelve.open('W_below_mxl_divH_ora_ab_ORAmesh.dat')
out['W_ab'] = wdTdz_ab_ab_weighted_b
out.close()
