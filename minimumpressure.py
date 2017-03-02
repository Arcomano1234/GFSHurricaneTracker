import numpy as np
import pygrib
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.ndimage as ndimage
from scipy.ndimage.filters import minimum_filter, maximum_filter
import os, time, datetime
from multiprocessing.dummy import Pool as ThreadPool 

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

test = '/data/studentdata/dkleist/tc/prctls15/'
color = '#ccffff'

def extractminpressure():
    path = '/data/studentdata/dkleist/tc/prctls15/' 
    directories = os.listdir(path)
    directories.sort() 
    for directory in directories:
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
                   x,y = m(lons,lats)
                   local_min, local_max = extrema(mslp, window=100)
                   ylows = y[local_min]
                   xlows = [local_min]     
                   lowvals = mslp[local_min] 
                   yoffset = 100000
                   dmin = yoffset
                   xyplotted = []
                   for hlx,hly,p in zip(xlows, ylows, lowvals):
                       if hlx < m.xmax and hlx > m.xmin and hly < m.ymax and hly > m.ymin:
                           dist = [np.sqrt((hlx-x0)**2+(hly-y0)**2) for x0,y0 in xyplotted]
                           if not dist or min(dist) > dmin:
                               if p < 1000 and p > 990:
                                  plt.plot(hlx,hly,'co')
                               elif p > 1000 and p < 1010:
                                  plt.plot(hlx,hly,'o',color = color)
                               elif p < 990 and p > 980:
                                  plt.plot(hlx,hly,'bo')
                               elif p < 980 and p > 970:
                                  plt.plot(hlx,hly,'go')
                               elif p < 970 and p > 960:
                                  plt.plot(hlx,hly,'yo')
                               elif p < 960 and p > 950:
                                  plt.plot(hlx,hly,'mo')
                               elif p < 950 and p > 940: 
                                  plt.plot(hlx,hly,'ro')
                               elif p < 940 and p > 930:
                                  plt.plot(hlx,hly,'ko')
                               elif p < 930: 
                                  plt.plot(hlx,hly,'wo')
                               else:
                                  pass
                               #plt.text(hlx+yoffset,hly+yoffset,repr(int(p)),fontsize=5,
                               #    ha='center',va='top',color='r',
                               #    bbox = bbox)
                               xyplotted.append((hlx,hly))
           #print p
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
           ax1 = figure.add_axes([0.14, 0.15, 0.2, 0.05])
           cmap = mpl.colors.ListedColormap([color,'c','b','g','y','m','r','k'])
           cmap.set_under('w')
           bounds = [1010,1000,990,980,970,960,950,940,930]
           norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
           cb = mpl.colorbar.ColorbarBase(ax1,cmap=cmap, norm=norm, boundaries = [880] + bounds, ticks = bounds,extend = 'min', orientation = 'Vertical')
           cb.set_label('hPa')
           plt.title(directory+'MSLP')

                   #plt.show()
           plotdir = '/homes/metogra/tarcomano/TropicalCycloneGFSJoaquin/python/'
           plt.savefig(plotdir+directory+'lowpressure.png')
           plt.close('all')

extractminpressure()
