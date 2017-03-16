# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:57:06 2017

@author: Tammas Loughran

Formula for streamfunction
              _p
psi = 2*pi*a | u_D* dp/g
            -0
             _p
   ~= 2*pi*a > u_D deltap/g
             ~0
"""
import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from windspharm.standard import VectorWind

# Equatorial Walkter circulation
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CMCC-CM','CNRM-CM5',
          'CSIRO-Mk3-6-0','GFDL-CM3','HadGEM2-A','IPSL-CM5A-LR','MIROC5',
          'MPI-ESM-MR','MRI-CGCM3','NorESM1-M','inmcm4']
elnino_years = [1982, 1987, 1991, 1997, 2002]
lanina_years = [1984, 1988, 1998, 1999, 2007]
cwd = os.getcwd()
amipdir = '/srv/ccrc/data48/z5032520/amip/'
a = 6.3710088E6 # meters
deltaphi = 2.0*np.pi #np.deg2rad(10)
g = 9.80665  # m/s^2

def locate3_in_list(values,alist):
    """Locates the indices of 'alist' where the first 3 instances of values
    in 'values' occurrs. These values should be Dec, Jan, Feb.
    
    Input
    values - a list of values to locate
    alist - the list to search
    
    Returns
    index - an array of indices
    """
    index = np.empty(0,dtype=int)
    for i in values:
        index = np.concatenate((index, np.where(alist==i)[0][-7:]), axis=0)
    return index

def plot_streamfunction(sfunc, ltude, pr, fname):
    fig = plt.figure(figsize=(8,4))
    levels = np.arange(-90.0E10,91E10,10E10)
    fill = plt.contourf(ltude,pr,sfunc,levels=levels,cmap='Spectral')
    plt.contour(ltude,pr,sfunc,levels=levels,colors='k')
    plt.plot([10,40],[980,980],'k',linewidth=5)
    plt.plot([280,320],[980,980],'k',linewidth=5)
    plt.plot([115,120],[980,980],'k',linewidth=5)
    cbar = plt.colorbar(fill)
    cbar.set_label('kg s$^{-1}$')
    
    ax = plt.gca()
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_yticks(np.array((1000,800,600,400,200,100)))
    ax.set_yticklabels(['1000','800','600','400','200','100'])
    ax.set_ylim([1000,100])
    ax.set_xticks(np.array((0,90,180,270,360)))
    ax.set_xticklabels(['0','90','180','270','360'])
    plt.ylabel('Pressure (hPa)')
    plt.xlabel('Longitude ($^{\circ}$E)')
    plt.savefig(fname+'.eps',format='eps')
    #plt.show()
    plt.close()


for model in models:
    print model
    os.chdir(cwd)
    os.chdir(model)
    # Identify ensemble members
    uensembles = glob.glob(amipdir+model+'/ua/*')
    vensembles = glob.glob(amipdir+model+'/va/*')
    n = 1
    for uens,vens in zip(uensembles,vensembles):
        # Load data
        ufiles = glob.glob(uens+'/ua_*')
        vfiles = glob.glob(vens+'/va_*')
        uanc = nc.MFDataset(ufiles,'r')
        vanc = nc.MFDataset(vfiles,'r')
        time = uanc.variables['time']
        time = nc.MFTime(time)
        time = nc.num2date(time[:],time.units)
        years = np.array([time[i].year for i in xrange(len(time))])
        months = np.array([time[i].month for i in xrange(len(time))])
        smr = (months==1)|(months==2)|(months==3)
        years[smr] = years[smr]-1
        lats = uanc.variables['lat'][:]
        lons = uanc.variables['lon'][:]
        p = uanc.variables['plev'][:]
        deltap = np.abs(np.diff(np.concatenate((p,[0])),n=1))
        ushape = uanc.variables['ua'].shape
        equator = (lats<5.)&(lats>-5)
        psi = np.ones((ushape[1],ushape[3]))*np.nan
        integral = np.zeros(ushape[3])
        for z in xrange(ushape[1]-1,-1,-1):
            u = uanc.variables['ua'][:,z,...]
            u[u>1000.] = 0
            u = np.array(u)
            v = vanc.variables['va'][:,z,...]
            v[v>1000.] = 0
            v = np.array(v)
            # Select and average El Nino years
            nino_index = locate3_in_list(elnino_years, years)
            nina_index = locate3_in_list(lanina_years, years)
            u_nino = u[nino_index].mean(axis=0)
            v_nino = v[nino_index].mean(axis=0)
            # Calculate the divergent wind
            uvwinds = VectorWind(u_nino,v_nino)
            u_D, _ = uvwinds.divergentcomponent()
            # Average over 5S to 5N
            u_D = u_D[equator,:].mean(axis=0)
            # Calculate streamfunction
            integral += u_D*deltap[z]/g
            psi[z,:] = a*deltaphi*integral
        
        # Plot the streamfunction
        pres = p/100.
        plot_streamfunction(psi, lons, pres, model+'_elninopsi_'+str(n))
        n += 1