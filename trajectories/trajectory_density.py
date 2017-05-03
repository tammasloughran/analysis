# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:12:39 2017

@author: Tammas Loughran
"""
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from subprocess import call
import glob
import scipy.stats as stats

def latlon2area(lat1,lat2,lon1,lon2):
    """Caculate the area of a latitude-longitude box on earth in km^2.
    """
    dlat = abs(np.sin(np.deg2rad(lat1)) - np.sin(np.deg2rad(lat2)))
    dlon = abs(lon1 - lon2)
    return (np.pi*6371.0088**2.)*dlat*dlon/180.

cwd = os.getcwd()
phase = 'lanina'
ntimes = 241
# Load the trajectories
def get_traj(cphase):
    trajdir = '/home/nfs/z5032520/traj3d/work/'
    os.chdir(trajdir+cphase)
    trajfiles = glob.glob('?????_heatwave_trajectories.nc')
    lons = np.ones((1,ntimes))*np.nan
    lats = np.ones((1,ntimes))*np.nan
    levs = np.ones((1,ntimes))*np.nan
    temp = np.ones((1,ntimes))*np.nan
    q = np.ones((1,ntimes))*np.nan
    for ifile in trajfiles:
        ncfile = nc.Dataset(ifile,'r')
        lons = np.concatenate((lons, np.squeeze(ncfile.variables['lon'][...])),axis=0)
        lats = np.concatenate((lats, np.squeeze(ncfile.variables['lat'][...])),axis=0)
        levs = np.concatenate((levs, np.squeeze(ncfile.variables['lev'][...])),axis=0)
        temp = np.concatenate((temp, np.squeeze(ncfile.variables['TEMP'][...])),axis=0)
        q = np.concatenate((q, np.squeeze(ncfile.variables['Q'][...])),axis=0)
    return lons, lats, levs, temp, q

lons, lats, levs, temp, q = get_traj('elnino')
alons, alats, alevs, atemp, aq = get_traj('lanina')


# Calculate potential temperature
Rcp = 0.286
theta = temp*(100000./levs)**Rcp
atheta = atemp*(100000./alevs)**Rcp
dxy = 2.
dlats = np.arange(-70.,1.,dxy)
dlons = np.arange(40.,201.,dxy)
# Create an array for the parcel density
def create_density(ibox):
    density = np.zeros((ntimes,len(dlats),len(dlons)))
    for t in xrange(ntimes):
        tlons = lons[ibox,t]
        tlats = lats[ibox,t]
        for y,dlat in enumerate(dlats):
            for x,dlon in enumerate(dlons):
                density[t,y,x] = ((tlons>=dlon)&
                                  (tlons<=dlon+dxy)&
                                  (tlats>=dlat)&
                                  (tlats<=dlat+dxy)).sum()
    prob = density/float(ibox.sum())*100
    return prob #np.ma.array(prob, mask=prob==0)


def create_track_density(ibox):
    density = np.ones((len(dlats),len(dlons)))*np.nan
    lont = lons[ibox,:]
    latt = lats[ibox,:]
    for y,dlat in enumerate(dlats):
        for x,dlon in enumerate(dlons):
            density[y,x] = ((lont>=dlon).any(axis=1)&
                            (lont<=dlon+dxy).any(axis=1)&
                            (latt>=dlat).any(axis=1)&
                            (latt<=dlat+dxy).any(axis=1)).sum()
    return np.ma.array(density,mask=density==0)


# Identify the trajectories that start within a given region.
def make_ibox(ur_box,size, lons, lats):
    cnrs = (ur_box[0],ur_box[0]+size,ur_box[1],ur_box[1]-size)
    ibox = (lons[:,0]>=cnrs[0])&(lons[:,0]<=cnrs[1])&(lats[:,0]<=cnrs[2])&(lats[:,0]>=cnrs[3])
    # Remove stationary trajectories that get stuck
    delta = np.sqrt((lons[:,0]-lons[:,-1])**2 + (lats[:,0]-lats[:,-1])**2)
    ibox[delta<10] = False
    return ibox


# Plot theta
def plot_theta(thetas, ntraj=0, outname='show'):
    fig = plt.figure()
    axes = fig.gca()
    bp = plt.boxplot(thetas, whis=0, showfliers=False)
    axes.set_xticks(np.arange(4,41,4))
    axes.set_xticklabels(np.arange(1,11,1))
    axes.set_ylim([0+273.15,60+273.15])
    axes.invert_xaxis()
    plt.xlabel('Days before heatwave')
    plt.ylabel(r'$\theta$ ($^{\circ}$C)')
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black',linestyle='-')
    #axes.spines['right'].set_visible(False)
    #axes.yaxis.set_ticks_position('left')
    #axes.spines['top'].set_visible(False)
    #axes.xaxis.set_ticks_position('bottom')
    plt.title('n='+str(ntraj))
    if outname=='show':
        plt.show()
    else:
        plt.savefig(outname, format='eps')
    plt.close()


def plot_levs(levls, ntraj=0, outname='show'):
    fig = plt.figure()
    axes = fig.gca()
    bp = plt.boxplot(levls, whis=0, showfliers=False)
    axes.set_xticks(np.arange(4,41,4))
    axes.set_xticklabels(np.arange(1,11,1))
    plt.xlabel('Days before heatwave')
    plt.ylabel(r'(hPa)')
    axes.invert_yaxis()
    axes.invert_xaxis()
    axes.set_yscale('log')
    axes.set_yticks(np.array((1000,800,600,400)))
    axes.set_yticklabels(['1000','800','600','400','200'])
    axes.set_ylim([1000,400])
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black',linestyle='-')
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.spines['top'].set_visible(False)
    axes.xaxis.set_ticks_position('bottom')
    plt.title('n='+str(ntraj))
    if outname=='show':
        plt.show()
    else:
        plt.savefig(outname, format='eps')
    plt.close()


def plot_density(ddata, ntraj=0, outname='show', dots=False, dot_lons=False, dot_lats=False):
    fig = plt.figure()
    xx,yy = np.meshgrid(dlons+1,dlats+1)
    m = Basemap(projection='mill',
                llcrnrlon=80.,llcrnrlat=-60.,
                urcrnrlon=200.,urcrnrlat=0.,
                resolution='l')
    xxx,yyy=m(xx,yy)
    #shade = m.pcolormesh(xxx,yyy,ddata,vmin=dens_levs[0],vmax=dens_levs[-1],cmap='viridis_r')
    #m.colorbar(shade)
    cont = m.contourf(xxx,yyy,ddata,
                      levels=dens_levs,
                      cmap='viridis_r',
                      #extend='max',
                      corner_mask=False)
    levels = [-1, 0.1]
    plt.contourf(xxx, yyy, ddata, levels=levels, colors='w')
    m.colorbar(cont)
    if dots:
        x, y = m(dot_lons[ibox],dot_lats[ibox])
        m.scatter(x,y,marker='.')
    m.drawcoastlines()
    plt.title('n='+str(ntraj))
    if outname=='show':
        plt.show()
    else:
        plt.savefig(outname, format='eps')
    plt.close()


def plot_density_wide(ddata, ntraj=0, outname='show', dots=False, dot_lons=False, dot_lats=False):
    fig = plt.figure()
    xx,yy = np.meshgrid(dlons,dlats)
    m = Basemap(projection='mill',
                llcrnrlon=40.,llcrnrlat=-70.,
                urcrnrlon=200.,urcrnrlat=0.,
                resolution='l')
    xxx,yyy=m(xx,yy)
    #shade = m.pcolormesh(xxx,yyy,ddata,vmin=dens_levs[0],vmax=dens_levs[-1],cmap='viridis_r')
    #m.colorbar(shade)
    cont = m.contourf(xxx,yyy,ddata,
                      levels=dens_levs,
                      cmap='viridis_r',
                      #extend='max',
                      corner_mask=False)
    levels = [-1, 0.1]
    plt.contourf(xxx, yyy, ddata, levels=levels, colors='w')
    m.colorbar(cont)
    if dots:
        x, y = m(dot_lons[ibox],dot_lats[ibox])
        m.scatter(x,y,marker='.')
    m.drawcoastlines()
    plt.title('n='+str(ntraj))
    if outname=='show':
        plt.show()
    else:
        plt.savefig(outname, format='eps')
    plt.close()


def random_traj(ibox1, ibox2):
    """Selct a random sample of 50 trajectories from within ibox1 and 
    ibox2 from elnino and lanina groups repectively.
    """
    nrtraj = 50
    rlons = np.ones((nrtraj,ntimes))*np.nan
    rlats = np.ones((nrtraj,ntimes))*np.nan
    rlevs = np.ones((nrtraj,ntimes))*np.nan
    rtemp = np.ones((nrtraj,ntimes))*np.nan
    rq = np.ones((nrtraj,ntimes))*np.nan
    for i in nrtraj:
        if (i%2==0):
            rn = np.random.randint(ibox1.sum())
            ri = np.where(ibox1==True)[0][rn]
            rlons[i,...] = lons[ri,...]
            rlats[i,...] = lats[ri,...]
            rlevs[i,...] = levs[ri,...]
            rtemp[i,...] = temp[ri,...]
            rq[i,...] = q[rn,...]
        else:
            rn = np.random.randint(ibox2.sum())
            ri = np.where(ibox2==True)[0][rn]
            rlons[i,...] = alons[ri,...]
            rlats[i,...] = alats[ri,...]
            rlevs[i,...] = alevs[ri,...]
            rtemp[i,...] = atemp[ri,...]
            rq[i,...] = aq[ri,...]
    return rlons, rlats, rlevs, rtemp, rq
    

def mannwhitneyu_2d(ad,bd):
    """Does a Mann-Whitney U test for 3 dimensional data allong the second axis.
    ad and bd are the two groups to test. 
    Returns us, the U statistic and pv the p-value for a two sided hypothesis test.
    """
    us = np.ones(ad.shape[1])*np.nan
    pv = us.copy()
    for y in xrange(ad.shape[1]):
        us[y], pv[y] = stats.mannwhitneyu(ad[:,y],bd[:,y], alternative='two-sided')
    return us, pv
    

def bootstrap_medians(ad,nsamp=10000):
    sample = np.ones((ad.shape))*np.nan
    medians = np.ones((nsamp,ad.shape[1]))*np.nan
    for n in xrange(nsamp):
        for i in xrange(ad.shape[0]):
            sample[i,...] = ad[np.random.randint(ad.shape[0]),...]
        medians[n,...] = np.median(sample, axis=0)
    lci = np.percentile(medians, 5, axis=0)
    uci = np.percentile(medians, 95, axis=0)
    return lci, uci
    
# Run
os.chdir(cwd)

size = 7.5
rnames = ['North','Northeast','Central-East','Southeast','Southwest']
regions = [(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)]
for i,ur_box in enumerate(regions):
    ibox = make_ibox(ur_box,size, lons, lats)
    num = ibox.sum()
    plot_theta(theta[ibox,2::6], 
               ntraj=num, 
               outname=rnames[i]+'_theta_'+phase+'.eps')
    plot_levs(levs[ibox,2::6]/100., 
              ntraj=num, 
              outname=rnames[i]+'levs'+phase+'.eps')
    prob = create_density(ibox)
    dens_levs = np.arange(0,11,1)
    for hour in [1,24,48,72]:
        dens_levs = np.arange(0,prob[hour,...].max(),1)
        plot_density(prob[hour,...], 
                     ntraj=num, 
                     outname=rnames[i]+'_'+phase+'_'+str(hour)+'.eps', 
                     dots=False)#, dot_lons=lons[:,hour], dot_lats=lats[:,hour])
    track = create_track_density(ibox)
    dens_levs = np.arange(0,track.max()+10,10)
    plot_density_wide(track, 
                 ntraj=num, 
                 outname='track_'+rnames[i]+'_'+phase+'.eps', 
                 dots=False)
                 
for i,ur_box in enumerate(regions):
    print rnames[i]
    ibox = make_ibox(ur_box,size, lons, lats)
    ibox2 = make_ibox(ur_box,size, alons, alats)
    num = ibox.sum()
    anum = ibox2.sum() 
    #rnlons, rnlats, rnlevs, rntemp, rnq = random_traj(ibox, ibox2)
    #_, thsig = mannwhitneyu_2d(theta[ibox,2::6],atheta[ibox2,2::6])
    #_, levsig = mannwhitneyu_2d(levs[ibox,2::6],alevs[ibox2,2::6])
    #print thsig<0.05
    #print levsig<0.05
    median = np.median(levs[ibox,2::6],axis=0)/100.
    amedian = np.median(alevs[ibox2,2::6],axis=0)/100.
    ninolci,ninouci = bootstrap_medians(levs[ibox,2::6])
    ninalci,ninauci = bootstrap_medians(alevs[ibox2,2::6])
    fig = plt.figure()
    axes = fig.gca()
    h1, = plt.plot(np.arange(0,40,1), median,color='r', label='El Nino n='+str(num))
    h2, = plt.plot(np.arange(0,40,1), amedian,color='b', label='La Nina n='+str(anum))
    plt.fill_between(np.arange(0,40,1), ninouci/100., ninolci/100., facecolor='coral',linewidth=0.0,alpha=1)
    plt.fill_between(np.arange(0,40,1), ninauci/100., ninalci/100., facecolor='cyan',linewidth=0.0,alpha=0.3)
    plt.legend(handles=[h1, h2])
    axes.set_xticks(np.arange(4,41,4))
    axes.set_xticklabels(np.arange(1,11,1))
    plt.xlabel('Days before heatwave')
    plt.ylabel(r'(hPa)')
    axes.invert_yaxis()
    axes.invert_xaxis()
    axes.set_yscale('log')
    axes.set_yticks(np.array((1000,900,800,600)))
    axes.set_yticklabels(['1000','900','800','600','400','200'])
    axes.set_ylim([1000,600])
    #plt.show()
    plt.savefig(rnames[i]+'_levels.svg',format='svg')
    plt.close()

# Plot all the trajectory points for animation
#for i in xrange(lons.shape[1]):
#    plot_density(density[i,...], outname=str(i).zfill(3), dot_lons=lons[:,i], dot_lats=lats[:,i])
#call(['ffmpeg', '-i %03d.png output.mov'])