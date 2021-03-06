import numpy as np
import scipy.ndimage
from scipy.ndimage.filters import minimum_filter, maximum_filter
import os, time, datetime
import matplotlib.pyplot as plt
import pygrib
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import multiprocessing as mp
import sys

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


path = '/data/studentdata/dkleist/tc/prctls15/2015093000/pgbq06.gfs.2015093000.grib2'
grbs = pygrib.open(path)

grb = grbs.select(name='MSLP (Eta model reduction)')[0]
mslp,lats,lons = grb.data(lat1=15,lat2=40,lon1=275,lon2=302)
mslp = mslp/100.
x, y = m(lons, lats)
local_min, local_max = extrema(mslp, window=75)
xlows = x[local_min]
ylows = y[local_min]
centerlat = lats[local_min]
centerlon = lons[local_min]
lowvals = mslp[local_min]

for hlx,hly in zip(centerlon,centerlat):
     lowerlon = hlx - 5
     upperlon = hlx + 5
     vort = []
     for grb in grbs:
         if grb.parameterName == 'Absolute vorticity' and grb.typeOfLevel == 'isobaricInhPa' and grb.level >= 70 and grb.level <= 1000:
             vorts,lats,lons=grb.data(lat1=hly,lat2=hly,lon1=lowerlon,lon2=upperlon)
             lons = np.squeeze(lons)
             vort.append(vorts)

     vort = np.squeeze(vort)      
     heights = np.array([[70,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000]])
     heights = np.squeeze(heights)
     fig = plt.figure()
     ax1 = fig.add_subplot(111)
     dx = 60 * np.cos(np.deg2rad(hly))
     dx = np.divide(dx,4)
     NormalizedVort = vort * 100000
     colortable =      [
                       "#ffffff",
                       "#ffff00",
                       "#ffee00",
                       "#ffdd00",
                       "#ffcc00",
                       "#ffbb00",
                       "#ffaa00",
                       "#ff9900",
                       "#ff8800",
                       "#ff7700",
                       "#ff6600",
                       "#ff5500",
                       "#ff4400",
                       "#ff3300",
                       "#ff2200",
                       "#ff1100",
                       "#ff0000",
                       ]
     vort_colormap = mpl.colors.ListedColormap(colortable)
     levels2 = np.arange(0,51,3)
     norm = mpl.colors.BoundaryNorm(levels2,19)
     ax1.contourf(lons,heights,NormalizedVort,levels=levels2,cmap=vort_colormap,norm=norm, extend='both')
     labelFontSize = "small"
     ax1.set_ylabel("Pressure [hPa]")
     ax1.set_yscale('log')
     ax1.set_ylim(10.*np.ceil(heights.max()/10.), heights.min()) # avoid truncation of 1000 hPa
     subs = [1,2,5]
     if heights.max()/heights.min() < 30.:
         subs = [1,2,3,4,5,6,7,8,9]
     y1loc = mpl.ticker.LogLocator(base=10., subs=subs)
     ax1.yaxis.set_major_locator(y1loc)
     fmt = mpl.ticker.FormatStrFormatter("%g")
     ax1.yaxis.set_major_formatter(fmt)
     for t in ax1.get_yticklabels():
         t.set_fontsize(labelFontSize)

     plt.show()


'''
grb = grbs.select(name='MSLP (Eta model reduction)')[0]
mslp,lats,lons = grb.data(lat1=15,lat2=40,lon1=275,lon2=302)
mslp = mslp/100.
x, y = m(lons, lats)
local_min, local_max = extrema(mslp, window=75)
xlows = x[local_min]
ylows = y[local_min]
centerlat = lats[local_min]
centerlon = lons[local_min]
lowvals = mslp[local_min]

        
grb = grbs.select(name='MSLP (Eta model reduction)')[0]
mslp,lats,lons = grb.data(lat1=15,lat2=40,lon1=275,lon2=302)
mslp = mslp/100.
x, y = m(lons, lats)
local_min, local_max = extrema(mslp, window=75)
xlows = x[local_min]
ylows = y[local_min]
centerlat = lats[local_min]
centerlon = lons[local_min]
lowvals = mslp[local_min]
print centerlon

heights = np.array([[70,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600]])

for hlx,hly in zip(centerlon,centerlat):
     lowerlon = hlx - 5
     upperlon = hlx + 5
     print upperlon
     print lowerlon
     fig = plt.figure()
     grbs.seek(0)
     selected_grbs = grbs(name='Absolute vorticity',typeOfLevel='isobaricInhPa',level= [70,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600])[0]
     ax1 = fig.add_subplot(111)
     absolutevort,lats,lons = selected_grbs.data(lat1=hly,lat2=hly,lon1=lowerlon,lon2=upperlon)
     NormalizedVort = absolutevort * 100000
     colortable = [
                  "#ffffff",
                  "#ffff00",
                  "#ffee00",
                  "#ffdd00",
                  "#ffcc00",
                  "#ffbb00",
                  "#ffaa00",
                  "#ff9900",
                  "#ff8800",
                  "#ff7700",
                  "#ff6600",
                  "#ff5500",
                  "#ff4400",
                  "#ff3300",
                  "#ff2200",
                  "#ff1100",
                  "#ff0000",
                  ]
     vort_colormap = mpl.colors.ListedColormap(colortable)
     levels2 = np.arange(0,51,30)
     print lons
     #lons,heights = np.meshgrid(lons,heights)
     print np.shape(NormalizedVort)
     print heights
     norm = mpl.colors.BoundaryNorm(levels2,19)
     ax1.contourf(lons,heights,NormalizedVort,levels=levels2,cmap=vort_colormap,norm=norm, extend='both') 
     ax1.set_yscale('log')
      
     plt.show()

minlowval = lowvals.min()
print minlowval
index = lowvals.index(minlowval)
print index
centerlat2 = lats[index]

def plotZM(data, x, y, plotOpt=None, modelLevels=None, surfacePressure=None):
    """Create a zonal mean contour plot of one variable
    plotOpt is a dictionary with plotting options:
      'scale_factor': multiply values with this factor before plotting
      'units': a units label for the colorbar
      'levels': use list of values as contour intervals
      'title': a title for the plot
    modelLevels: a list of pressure values indicating the model vertical resolution. If present,
        a small side panel will be drawn with lines for each model level
    surfacePressure: a list (dimension len(x)) of surface pressure values. If present, these will
        be used to mask out regions below the surface
    """
    # explanation of axes:
    #   ax1: primary coordinate system latitude vs. pressure (left ticks on y axis)
    #   ax2: twinned axes for altitude coordinates on right y axis
    #   axm: small side panel with shared y axis from ax2 for display of model levels
    # right y ticks and y label will be drawn on axr if modelLevels are given, else on ax2
    #   axr: pointer to "right axis", either ax2 or axm

    if plotOpt is None: plotOpt = {}
    labelFontSize = "small"
    # create figure and axes
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # scale data if requested
    scale_factor = plotOpt.get('scale_factor', 1.0)
    pdata = data * scale_factor
    # determine contour levels to be used; default: linear spacing, 20 levels
    clevs = plotOpt.get('levels', np.linspace(data.min(), data.max(), 20))
    # map contour values to colors
    norm=matplotlib.colors.BoundaryNorm(clevs, ncolors=256, clip=False)
    # draw the (filled) contours
    contour = ax1.contourf(x, y, pdata, levels=clevs, norm=norm) 
    # mask out surface pressure if given
    if not surfacePressure is None: 
        ax1.fill_between(x, surfacePressure, surfacePressure.max(), color="white")    
    # add a title
    title = plotOpt.get('title', 'Vertical cross section')
    ax1.set_title(title)
    # add colorbar
    # Note: use of the ticks keyword forces colorbar to draw all labels
    fmt = matplotlib.ticker.FormatStrFormatter("%g")
    cbar = fig.colorbar(contour, ax=ax1, orientation='horizontal', shrink=0.8,
                        ticks=clevs, format=fmt)
    cbar.set_label(plotOpt.get('units', ''))
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(labelFontSize)
    # set up y axes: log pressure labels on the left y axis, altitude labels
    # according to model levels on the right y axis
    ax1.set_ylabel("Pressure [hPa]")
    ax1.set_yscale('log')
    ax1.set_ylim(10.*np.ceil(y.max()/10.), y.min()) # avoid truncation of 1000 hPa
    subs = [1,2,5]
    if y.max()/y.min() < 30.:
        subs = [1,2,3,4,5,6,7,8,9]
    y1loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
    ax1.yaxis.set_major_locator(y1loc)
    fmt = matplotlib.ticker.FormatStrFormatter("%g")
    ax1.yaxis.set_major_formatter(fmt)
    for t in ax1.get_yticklabels():
        t.set_fontsize(labelFontSize)
    # calculate altitudes from pressure values (use fixed scale height)
    z0 = 8.400    # scale height for pressure_to_altitude conversion [km]
    altitude = z0 * np.log(1015.23/y)
    # add second y axis for altitude scale 
    ax2 = ax1.twinx()
    # change values and font size of x labels
    ax1.set_xlabel('Latitude [degrees]')
    xloc = matplotlib.ticker.FixedLocator(np.arange(-90.,91.,30.))
    ax1.xaxis.set_major_locator(xloc)
    for t in ax1.get_xticklabels():
        t.set_fontsize(labelFontSize)
    # draw horizontal lines to the right to indicate model levels
    if not modelLevels is None:
        pos = ax1.get_position()
        axm = fig.add_axes([pos.x1,pos.y0,0.02,pos.height], sharey=ax2)
        axm.set_xlim(0., 1.)
        axm.xaxis.set_visible(False)
        modelLev = axm.hlines(altitude, 0., 1., color='0.5')
        axr = axm     # specify y axis for right tick marks and labels
        # turn off tick labels of ax2
        for t in ax2.get_yticklabels():
            t.set_visible(False)
        label_xcoor = 3.7
    else:
        axr = ax2
        label_xcoor = 1.05
    axr.set_ylabel("Altitude [km]")
    axr.yaxis.set_label_coords(label_xcoor, 0.5)
    axr.set_ylim(altitude.min(), altitude.max())
    yrloc = matplotlib.ticker.MaxNLocator(steps=[1,2,5,10])
    axr.yaxis.set_major_locator(yrloc)
    axr.yaxis.tick_right()
    for t in axr.yaxis.get_majorticklines():
        t.set_visible(False)
    for t in axr.get_yticklabels():
        t.set_fontsize(labelFontSize)
    # show plot
    plt.show()
'''
