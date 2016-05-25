#!/usr/bin/env python

from netCDF4 import Dataset
import plot_tool as pf
import numpy as np
from numpy import inf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import maskoceans
import os
import pickle
import shelve
import copy
import math

def shift_data(data,lon):
  """
  Shifts data on a map by 180 degree longitude.
  Input: data = your array and longitude array.
  """
  if len(data.shape) = 3:
  	data1 = data[:,:,:len(lon)/2+1]
  	data2 = data[:,:,len(lon)/2+1:]
  	data_180 = np.concatenate((data2,data1),axis=2)
  if len(data.shape) = 4:
  	data1 = data[:,:,:,:len(lon)/2+1]
  	data2 = data[:,:,:,len(lon)/2+1:]
  	data_180 = np.concatenate((data2,data1),axis=2)
  else:
  	print "I don't understand your data"
  return data_180

def shift_lon(lon):
  """
  Shifts you lon so that you have it to plot
  """
  lon1 = lon[:len(lon)/2+1]
  lon2 = lon[len(lon)/2+1:]
  for i in range(len(lon2)):
    lon2[i] = lon2[i]-360
  lon_180 = np.r_[lon2,lon1]
  return lon_180

def extract_months(data,month,ly):
  """
  Extracts relevant months from timeseries of your experiment
  Input: data = your array, month = which month do you want to extract
  (may, jun, jul, aug)
  ly = how many years are there?
  """
  ms = [31,30,31,31]
  mos1 = ['may','jun','jul','aug']
  days_year = sum( ms )
  for i in range( ly ):
    if i <= 0:
      dummy =  data[i*(days_year):sum(ms[:mos1.index(month)+1]),:,:]
      m_all = dummy
    else:
      dummy =  data[i*(days_year):i*days_year +
		    sum(ms[:mos1.index(month)+1]),:,:]
      m_all = np.concatenate((m_all[:,:,:],dummy[:,:,:]),axis=0)
  return m_all

def extract_years(data,lat,lon,leny,ny):
  """
  Extracts years from timeseries of your experiment
  Input: data = your array, leny = how long is a year (ts)
  ny = how many years are there?
  """
  years = np.zeros((ny,leny,len(lat),len(lon)))
  for i in range( ny ):
    years[i,:,:,:] =  data[i*leny:i*leny+leny,:,:]
  return years

def average_month(data,lm,ly):
  """
  Averages a daily time series with multiple months 
  over the month.
  Input: data = your array, lm = length of the month,
  ly = how many years?
  """
  months = np.zeros((ly,lm))
  for i in range(ly):
    months[i,:] = data[i*lm:i*lm+lm]
  ave_months = np.nanmean(months, axis=0)
  return ave_months

def find_range(lon,lat,lon_min,lon_max,lat_min,lat_max):
  """
  Finds the lat and lon indices for the box you want.
  Assums that lat goes from 90 to -90
  and lon from -180 to 180.
  Returns lomin, lomax, lamin, lamax.
  """
  lomin = np.min(np.where(lon > lon_min))
  lomax = np.max(np.where(lon < lon_max))
  lamax = np.max(np.where(lat > lat_min)) # of course this depends
  lamin = np.min(np.where(lat < lat_max)) # on the starting point
  return lomin, lomax + 1, lamin, lamax + 1

def calc_mean(data):
  """
  Calculate time mean of your data.
  """
  mean = np.nanmean(data[:,:,:],axis=0) #contains infinity 
  if np.isinf(mean).any():
    mean[mean == inf] = np.nan
  return mean

def sellonlat(lon,lat,lon_min,lon_max,lat_min,lat_max,data):
  """
  Finds the lat and lon indices for the box you want.
  Assums that lat goes from 90 to -90
  and lon from -180 to 180.
  Then selects the data for the box you want.
  """
  lomin = np.min(np.where(lon > lon_min))
  lomax = np.max(np.where(lon < lon_max))
  lamax = np.max(np.where(lat > lat_min))
  lamin = np.min(np.where(lat < lat_max))
  cut_data = data[:,lamin:lamax+1,lomin:lomax+1]
  return cut_data
