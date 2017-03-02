import numpy as np
import pygrib
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.ndimage as ndimage
from scipy.ndimage.filters import minimum_filter, maximum_filter
import os, time, datetime
import multiprocessing as mp 
import sys

# TODO get lats and lons of center to a text file to work with 
#f = open("log.testmultiprocessor.py",'w')
#sys.stdout = f

def extrema(mat,mode='wrap',window=10):
    """find the indices of local extrema (min and max)
    in the input array."""
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)

m=Basemap(projection='mill',lat_ts=10,llcrnrlon=270, \
  urcrnrlon=310,llcrnrlat=10,urcrnrlat=50, \
  resolution='i') 

bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.75))

color = '#ccffff'
orange = '#ff6600'
path = '/data/studentdata/dkleist/tc/prctls15/'

def hurricanetracker(directory):
    if directory != "cp.txt" and directory != "unpack.sh":
        figure = plt.figure(figsize=(12,8))
        modelrun = directory
        pathfiles = path+directory+'/'
        files = os.listdir(pathfiles)
        for file in files:
            if file[-5:] == "grib2":
                grib = path+directory+'/'+file
                print file
                grbs = pygrib.open(grib)
                grb = grbs.select(name='MSLP (Eta model reduction)')[0]
                mslp,lats,lons = grb.data(lat1=15,lat2=45,lon1=275,lon2=302)
                mslp = mslp/100.
                x, y = m(lons, lats)
                local_min, local_max = extrema(mslp, window=75)
                xlows = x[local_min]
                ylows = y[local_min]
                lowvals = mslp[local_min]
                centerlat = lats[local_min]
                centerlon = lons[local_min] 
               
                grbs.rewind()
                grb2 = grbs.select(name = 'Temperature', typeOfLevel ='isobaricInhPa', level = 850)[0]
                temp2,lats,lons = grb2.data(lat1=15,lat2=45,lon1=275,lon2=302)
                tempvals = temp2[local_min]
                
                grbs.rewind()
                grb3 = grbs.select(name='Absolute vorticity',typeOfLevel='isobaricInhPa',level=850)[0]
                vort850,lats,lons = grb3.data(lat1=15,lat2=45,lon1=275,lon2=302) 
                vortvals = vort850[local_min]
                yoffset = 100000
                dmin = yoffset
                xyplotted = []
                for hlx,hly,p,temp,vort in zip(xlows, ylows, lowvals, tempvals, vortvals):
                     cutoff = 0
                     if temp > 286:
                        cutoff = 1 + cutoff
                     elif temp > 283 and temp < 286:
                        cutoff = 0.5 + cutoff
                     else:
                        pass
                     if vort > 0.00015:
                        cutoff = 1 + cutoff
                     elif vort > 0.0001 and vort < 0.00015: 
                        cutoff = cutoff + 0.5
                     else:
                        pass
                     #print cutoff       
                     if hlx < m.xmax and hlx > m.xmin and hly < m.ymax and hly > m.ymin and cutoff >= 1:
                         dist = [np.sqrt((hlx-x0)**2+(hly-y0)**2) for x0,y0 in xyplotted]
                         if not dist or min(dist) > dmin:
                             if p < 1000 and p > 990:
                                plt.plot(hlx,hly,'co')
                             elif p > 1000 and p < 1020:
                                plt.plot(hlx,hly,'o',color = color)
                             elif p < 990 and p > 980:
                                plt.plot(hlx,hly,'bo')
                             elif p < 980 and p > 970:
                                plt.plot(hlx,hly,'go')
                             elif p < 970 and p > 960:
                                plt.plot(hlx,hly,'yo')
                             elif p < 960 and p > 950:
                                plt.plot(hlx,hly,'o',color = orange)
                             elif p < 950 and p > 940: 
                                plt.plot(hlx,hly,'ro')
                             elif p < 940 and p > 930:
                                plt.plot(hlx,hly,'mo')
                             elif p < 930: 
                                plt.plot(hlx,hly,'wo')
                             else:
                                pass
                             xyplotted.append((hlx,hly))
                     else:
                        pass
                     #TODO eventually put return xyplotted so I can get the lat and lon centers
        hlx = []
        hly = []
        p = [] 
        del hlx
        del hly
        del p  
        #cs = m.pcolormesh(x,y,mslp,shading='flat',cmap=plt.cm.jet)
        m.drawcoastlines()
        m.fillcontinents()
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,120.,5.),labels=[1,0,0,0])
        m.drawmeridians(np.arange(-180.,180.,5.),labels=[0,0,0,1])
        plt.title(directory+'MSLP With Relocation (Operational Version)')
        ax1 = figure.add_axes([0.8, 0.3, 0.02, 0.4])
        cmap = mpl.colors.ListedColormap(['m', 'r', orange, 'y', 'g', 'b', 'c'])
        cmap.set_under('w')
        cmap.set_over(color)
        bounds = [930, 940, 950, 960, 970, 980, 990, 1000]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        cb2 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,
                                        norm=norm,
                                   # to use 'extend', you must
                                 # specify two extra boundaries:
                                        boundaries=[900] + bounds + [1010],
                                        extend='both',
                                        ticks=bounds,  # optional
                                        spacing='uniform',
                                        orientation='vertical')
        cb2.set_label('hPa')
        plotdir = '/homes/metogra/tarcomano/TropicalCycloneGFSJoaquin/python/'
        plt.savefig(plotdir+directory+'lowpressure.png')
        plt.close('all')

def multip_handler():
    directories = os.listdir(path)
    directories.sort()
    p = mp.Pool(4)
    p.map(hurricanetracker,directories)
    p.close()
    p.join()

if __name__ == '__main__':
   multip_handler()

