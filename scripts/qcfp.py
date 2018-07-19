# standard modules
import datetime
import logging
import os
import zipfile
# 3rd party modules
import matplotlib.pyplot as plt
import numpy
import pandas
#from scipy.misc.pilutil import imread
# PFP modules
import constants as c
import qcgf
import qcio
import qcts
import qcutils

#rclim = sys.getrecursionlimit()
#print 'rclim = ', rclim
#sys.setrecursionlimit(10000)

# Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P., 2015:
# A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP),
# Geosci. Model Dev., 8, 3695-3713.
import calc_footprint_FFP_climatology as calcfootNK
# Kormann, R. and Meixner, F.X., 2001: An analytical footprint model for non-neutral stratification.
# Boundary-Layer Meteorology 99: 207. https://doi.org/10.1023/A:1018991015119 and Neftel, A., Spirig, C.,
# Ammann, C., 2008: Application and test of a simple tool for operational footprint evaluations. Environmental
# Pollution 152, 644-652.
import calc_footprint_FKM_climatology as calcfootKM

logger = logging.getLogger("pfp_log")

# constant for converting degrees to radiant
c_d2r = numpy.pi/180.0
# constant to convert the footprint area from x,y in lon,lat coordinates at the tower site
onedegree = 6378100.0 * c_d2r # distance in m for 1 deg of latitude

# Coordinate steps in footprint process
def footprint_main(cf, mode):
    """
    This script reads data from a PyFluxPro .nc file and processes the data for:
    (1) Kormann&Meixner uses input (single time)        zm,z0,ustar,umean,L,sigmav,wind_dir
    (2) Kljun et al. uses input (vectors for all times) zm,z0,ustar,umean,L,sigmav,wind_dir,Habl
        so Natascha's FFP also needs the height of the boundary layer ===> currently ERAI data got Habl,
        and ACCESS got Habl00 ... Habl22
    === > input for Kormann & Meixner and Natascha Kljun's footprint climatology
    (a) PyFluxPro L3 netcdf file
    (b) ERAI/ACCESS netcdf file
    === > output for the climatology
    (a) daily footprint climatology
    (b) monthly footprint climatology
    (c) annual footprint climatology
    (d) special time set in controlfile for footprint climatology
    GOAL: Footprint climatology can be done on a set time in controlfile
          calculating Kljun et al., 2015 and Kormann and Meixner, 2001 footprint
    DONE: set time in controlfile, special, daily, monthly and annual
          Kljun et al. (2015) footprint
          Kormann and Meixner (2001) footprint
          save footprint fields in netcdf file
    Still to do: calculate Habl if not exist
                 plot footprint in googleEarth kml format
    C.M.Ewenz, 10 June 2018
               21 June 2018 (corrections to monthly indexing)
               29 June 2018 (kml file, single time stamp)
    """
    logger.info(' Read input data files ...')
    # get the L3 data
    ds = get_footprint_data_in(cf, mode)
    ldt = ds.series["DateTime"]["Data"]
    # get the configuration data for the footprint
    d = get_footprint_cfg(cf, ds)
    logger.info(' Starting footprint calculation ...')
    list_StDate, list_EnDate = create_index_list(cf, d, ldt)
    logger.info(' Starting footprint climatology calculation ...')
    # !!! Prepare Output netcdf file !!!
    # Set initial x,y Variables for output
    xout = numpy.linspace(d["xmin"], d["xmax"], d["nx"] + 1)
    yout = numpy.linspace(d["ymin"], d["ymax"], d["nx"] + 1)
    lat0 = float(d["latitude"])
    lon0 = float(d["longitude"])
    lat = lat0 + yout / onedegree
    lon = lon0 + xout / (numpy.cos(lat0*c_d2r) * onedegree)

    # - Initialise output netcdf file and write x,y grid into file as xDistance and yDistance from the tower
    nc_name = d["file_out"]
    #print 'nc_name = ',nc_name
    nc_file = qcio.nc_open_write(nc_name)
    # create the x and y dimensions.
    nc_file.createDimension('longitude', len(lon))
    nc_file.createDimension('latitude', len(lat))
    # create time dimension (record, or unlimited dimension)
    nc_file.createDimension('time', None)
    # create number of footprints in climatology dimension (record, or unlimited dimension)
    nc_file.createDimension('dtime', None)
    nc_file.createDimension('num', None)
    # Define coordinate variables, which will hold the coordinate information, x and y distance from the tower location.
    X = nc_file.createVariable('longitude', "d", ('longitude',))
    Y = nc_file.createVariable('latitude', "d", ('latitude',))
    # Define time variable and number of footprints variable at each time
    tx = nc_file.createVariable('dtime', "d", ('dtime',))
    num = nc_file.createVariable('num', "d", ('num',))
    # Assign units attributes to coordinate var data, attaches text attribute to coordinate variables, containing units.
    X.units = 'degree'
    Y.units = 'degree'
    # write data to coordinate vars.
    X[:] = lon
    Y[:] = lat
    # create the sumphi variable
    phi = nc_file.createVariable('sumphi', "d", ('time', 'longitude', 'latitude'))
    # set the units attribute.
    phi.units = ' '
    # === General inputs for FFP
    zmt = d["zm_d"]
    domaint = [d["xmin"], d["xmax"], d["ymin"], d["ymax"]]
    nxt = d["nx"]
    rst = None #[90.] #None #[20.,40.,60.,80.]
    # if plotting to screen      is requested then iplot = 1
    # if plotting to googleEarth is requested then iplot = 2
    iplot = int(cf['General']['iplot'])

    #PlotWidth = float(cf['General']['PlotWidth'])
    #PlotHeight = float(cf['General']['PlotHeight'])

    # do we want to export images in kml format?
    if iplot == 2:  # kml - format header

        kmlname = d["site_name"] + '_' + mode + '_fp' + '.kml'
        kml_name_path = d["plot_path"] +kmlname
        fi = open(kml_name_path, 'w')
        fi.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fi.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        fi.write("<Folder>\n")
        fi.write("  <name>" + d["site_name"] + "</name>")
        # make sure GE zooms in
        fi.write('  <LookAt>\n')
        fi.write('    <longitude>'+str(d["longitude"])+'</longitude>\n')
        fi.write('    <latitude>'+str(d["latitude"])+'</latitude>\n')
        fi.write('    <altitude>'+str(d["footprint_size"])+'</altitude>\n')
        fi.write('    <range>'+str(d["footprint_size"])+'</range>\n')
        fi.write('    <tilt>0</tilt>\n')
        fi.write('    <heading>0</heading>\n')
        fi.write('    <altitudeMode>relativeToGround</altitudeMode>\n')
        fi.write('  </LookAt>\n')
        # ------------------------

    # After deciding which climatology is done, let's do it!

    irun = -1
    for i in range(0, len(list_StDate)):
        irun = irun+1
        # original CE code
        # ===========================================================================
        #umeant   = list(temp_df["Ws"][list_StDate[i]:list_EnDate[i]].astype(float))
        #olt      = list(temp_df["L"][list_StDate[i]:list_EnDate[i]].astype(float))
        #sigmavt  = list(temp_df["StdV"][list_StDate[i]:list_EnDate[i]].astype(float))
        #ustart   = list(temp_df["ustar"][list_StDate[i]:list_EnDate[i]].astype(float))
        #wind_dirt= list(temp_df["Wd"][list_StDate[i]:list_EnDate[i]].astype(float))
        #ldt      = list(temp_df["xlDateTime"][list_StDate[i]:list_EnDate[i]].astype(float))
        #z0t      = list(temp_df["z0"][list_StDate[i]:list_EnDate[i]].astype(float))
        # PRI suggested changes
        # 1) get the required series from the data structure as masked arrays
        # 2) get a mask that is True where the elements of 1 or more series are masked
        # 3) apply the comosite mask to the series
        # 4) compress the masked series to remove masked elements
        # get the start and end indices
        si = list_StDate[i]
        ei = list_EnDate[i]
        # get the series as masked arrays
        umeant, _, _ = qcutils.GetSeriesasMA(ds, "Ws", si=si, ei=ei)
        olt, _, _ = qcutils.GetSeriesasMA(ds, "L", si=si, ei=ei)
        sigmavt, _, _ = qcutils.GetSeriesasMA(ds, "V_Sd", si=si, ei=ei)
        ustart, _, _ = qcutils.GetSeriesasMA(ds, "ustar", si=si, ei=ei)
        wind_dirt, _, _ = qcutils.GetSeriesasMA(ds, "Wd", si=si, ei=ei)
        z0t, _, _ = qcutils.GetSeriesasMA(ds, "z0", si=si, ei=ei)
        ht, _, _ = qcutils.GetSeriesasMA(ds, "Habl", si=si, ei=ei)
        # get a composite mask over all variables
        mask_all = numpy.ma.getmaskarray(ustart)
        for item in [umeant, olt, sigmavt, wind_dirt, z0t, ht]:
            mask_item = numpy.ma.getmaskarray(item)
            mask_all = numpy.ma.mask_or(mask_all, mask_item)
        # and then apply the composite mask to all variables and remove masked elements
        umeant = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, umeant)))
        olt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, olt)))
        sigmavt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, sigmavt)))
        ustart = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ustart)))
        wind_dirt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, wind_dirt)))
        z0t = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, z0t)))
        ht = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ht)))
        if mode == "kljun":
            #ht       = list(temp_df["Habl"][list_StDate[i]:list_EnDate[i]].astype(float))
            # now run the footprint climatology for a day/month/year at a time ## calcfoot.
            FFP = calcfootNK.FFP_climatology (zm=zmt,z0=z0t,umean=umeant,h=ht,ol=olt,sigmav=sigmavt,ustar=ustart,\
                                              wind_dir=wind_dirt,domain=domaint,dx=None,dy=None,nx=nxt,ny=None,\
                                            rs=rst,rslayer=0,smooth_data=1,crop=False,pulse=None,verbosity=2)

            x              = FFP['x_2d']
            y              = FFP['y_2d']
            f              = FFP['fclim_2d']
            num[irun]      = FFP['n']

            #tx[irun] = int(ldt[0])
            phi[irun,:,:] = f

            fmax=numpy.amax(f)
            f=f/fmax

        elif mode == "kormei":
            #logger.info(' Kormann & Meixner footprint not yet implemented.')
            FKM = calcfootKM.FKM_climatology(zm=zmt, z0=z0t, umean=umeant, ol=olt, sigmav=sigmavt, ustar=ustart,\
                                             wind_dir=wind_dirt, domain=domaint, dx=None, dy=None, nx=nxt, ny=None, \
                     rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], rslayer=0,\
                     smooth_data=1, crop=False, pulse=None, verbosity=0)

            x              = FKM['x_2d']
            y              = FKM['y_2d']
            f              = FKM['fclim_2d']
            num[irun]      = FKM['n']

            #tx[irun] = int(ldt[0])
            phi[irun,:,:] = f

            fmax=numpy.amax(f)
            f=f/fmax
        else:
            msg = " Unrecognised footprint type " + str(mode)
            logger.error(msg)
            return
        # ====================================================================================================
        # get the default plot width and height
        clevs = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        imagename = qcutils.get_keyvaluefromcf(cf,["General"],'OzFlux_area_image')
        if iplot == 1:
            # get the image filename if plotting if OzFlux_area_image in cf['General'].keys():
            plotphifield(x, y, ldt[si], ldt[ei], imagename, d["site_name"], mode)
            #plot = plotphifield(temp_df["xlDateTime"][list_StDate[i]],x,y,\
                                #temp_df["xlDateTime"][list_StDate[i]],\
                                #temp_df["xlDateTime"][list_EnDate[i]],\
                                #f,imagename,d["site_name"],mode,clevs)
        elif iplot == 2:
            #pytz.timezone(d["timezone"]).localize(datetime.datetime(2011,1,1)).strftime('%z')
            # write a kml file with the image data
            write_kml(lon, lat, ldt[0], ldt[-1], f, d["site_name"], mode, clevs, fi, d["plot_path"])

    if iplot == 2:
        # write the footer of the kml file and close the file
        fi.write("</Folder>\n")
        fi.write('</kml>\n')
        fi.close()
        # create a kmz file out of the kml file
        #print
        #print '----------------'
        #print 'Zipping KMZ...'
        #print d["plot_path"]
        cwd = os.getcwd()
        #print cwd
        os.chdir(d["plot_path"])
        kmzname = kmlname.replace(".kml", ".kmz")
        msg = " Creating KMZ file " + kmzname
        logger.info(msg)
        #print kmlname,kmzname
        #print os.listdir('.')
        plotlist = [p for p in os.listdir('.') if p.endswith(".png")]
        #print plotlist
        compression = zipfile.ZIP_DEFLATED
        zf = zipfile.ZipFile(kmzname, mode='w')
        zf.write(kmlname, compress_type=compression)
        os.remove(kmlname)
        for f in plotlist:
            zf.write(f, compress_type=compression)
            os.remove(f)
        zf.close()
        os.chdir(cwd)
        #print
        #print 'Zipped:', kmzname
        #print

    # ================================================================
    msg = " Finished " + str(mode) + " footprint writing"
    logger.info(msg)
    msg = " Closing netcdf file " + str(nc_name)
    logger.info(msg)
    nc_file.close()
    # ================================================================

def write_kml(lon, lat, zt1, zt2, data, station, mode, clevs, fi, plot_path):

    plot_in='Footprint_'+ mode + zt1.strftime("%Y%m%d%H%M") +'.png'
    plotname=plot_path + plot_in

    width = 5
    height = width * data.shape[0]/data.shape[1]
    fig=plt.figure(figsize=(width,height))
    fig.add_axes([0,0,1,1])

    plt.contourf(data,clevs,cmap=plt.get_cmap('hsv'),alpha=0.5)

    plt.axis('off')
    plt.savefig(plotname,transparent=True)

    # get the lat/lon bounds of the area
    lon1 = lon[0]
    lon2 = lon[-1]
    lat1 = lat[0]
    lat2 = lat[-1]

    # Hopefully the file was opened properly and the header written
    fi.write('<GroundOverlay>\n')
    fi.write('  <name>'+station+str(zt2)+'</name>\n')
    fi.write('  <bgColor>8fffffff</bgColor>\n')
    fi.write('  <Icon>\n')
    fi.write('    <href>'+plot_in+'</href>\n')
    fi.write('  </Icon>\n')
    fi.write('  <TimeSpan>\n')
    fi.write('    <begin>'+str(zt1)+'</begin>\n')
    fi.write('    <end>'+str(zt2)+'</end>\n')
    fi.write('  </TimeSpan>\n')
    fi.write('  <altitude>0.0</altitude>\n')
    fi.write('  <altitudeMode>clampToGround</altitudeMode>\n')
    fi.write('  <LatLonBox>\n')
    fi.write('    <north>'+str(lat1)+'</north>\n')
    fi.write('    <south>'+str(lat2)+'</south>\n')
    fi.write('    <east>'+str(lon1)+'</east>\n')
    fi.write('    <west>'+str(lon2)+'</west>\n')
    fi.write('    <rotation>0.0</rotation>\n')
    fi.write('  </LatLonBox>\n')
    fi.write('</GroundOverlay>\n')

def plotphifield(x, y, zt1, zt2, imagename, station, add):

    text = 'Footprint ' + station + ' ' + zt1.strftime("%Y%m%d%H%M") + '  to  ' + zt2.strftime("%Y%m%d%H%M")
    plotname='plots/Footprint_'+ add + zt1.strftime("%Y%m%d%H%M") + '.jpg'
    x_ll = x[0,0]   #xllcorner #-250
    x_ur = x[-1,-1] #xurcorner # 250
    y_ll = y[0,0]   #yllcorner #-250
    y_ur = y[-1,-1] #yurcorner # 250
    # create figure and axes instances
    fig = plt.figure(figsize=(10,10))
    # contour levels
    plt.title(text)
    if imagename != 'None':
        #img = imread(imagename)
        img = plt.imread(imagename)
        plt.imshow(img, zorder=0, extent=[x_ll, x_ur, y_ll, y_ur])
    plt.savefig(plotname)
    fig.show()

def xldate_to_datetime(xldate):
    tempDate = datetime.datetime(1899,12,30)
    deltaDays = datetime.timedelta(days=int(xldate))
    secs = (int((xldate%1)*86400))
    deltaSeconds = datetime.timedelta(seconds=secs)
    TheTime = (tempDate + deltaDays + deltaSeconds )
    #return TheTime.strftime("%d-%m-%Y %H:%M")
    return TheTime.strftime("%Y-%m-%dT%H:%M")

def z0calc(zm,LM,U_meas,UStar):
    # aerodynamic roughness length
    # Psi functions according to Dyer (1974)
    # a) create positive and negative LM masks
    LMp = numpy.ma.masked_where(LM <  float(0),LM)
    LMn = numpy.ma.masked_where(LM >= float(0),LM)
    # Calculate z0 assuming logarithmic wind profile
    #          === functions are from Kormann and Meixner (2001) (Eqs. 31 to 35)
    #b) for stable conditions, linear
    FIp = 5.0 * zm/LMp
    # c) for unstable conditions
    zeta = (1.0-16.0*zm/LMn)**(0.25)
    FIn = -2.0*numpy.log(0.5*(1.0+zeta))-numpy.log(0.5*(1.0+zeta*zeta))+2.0*numpy.arctan(zeta)-0.5*c.Pi
    # d) put both parts together again
    #FI = numpy.ma.mask_or(FIp,FIn)
    # d1) fill positive and negative Fn masks
    FIp = numpy.ma.filled(FIp,float(0))
    FIn = numpy.ma.filled(FIn,float(0))
    FI  = FIp+FIn
    # e) determine
    alpha = U_meas * 0.4 / UStar - FI
    # f) finally calculate the roughness length
    ZNull = zm / numpy.exp(alpha)
    #!#            === functions derived from Leclerc and Foken 2015 book, page 61 after Hogstroem, 1988
    #!# b) for stable conditions, linear
    #!FIp = -6.0 * zm/LMp
    #!# c) for unstable conditions
    #!zeta = (1.0+19.3*zm/LMn)**0.25
    #!temp = 0.125*(1.0+zeta*zeta)*(1.0+zeta)*(1.0+zeta)
    #!FIn = numpy.log(temp)-2.0*numpy.arctan(zeta)+0.5*c.Pi
    #!# d) put both parts together again
    #!#FI = numpy.ma.mask_or(FIp,FIn,copy=True)
    #!# d1) fill positive and negative Fn masks
    #!FIp = numpy.ma.filled(FIp,float(0))
    #!FIn = numpy.ma.filled(FIn,float(0))
    #!FI  = FIp+FIn
    #!# e) determine
    #!alpha = U_meas * 0.4 / UStar + FI
    #!# f) finally calculate the roughness length
    #!ZNull = zm / numpy.exp(alpha)
    #!#            ===
    #set a lower limit for z0 to avoid numeric problems
    ZNull = numpy.ma.masked_where(ZNull<0.0001,ZNull)
    ZNull = numpy.ma.filled(ZNull,0.0001)
    return ZNull

#def BLH(ol, ustar, dt, lat, zm):
def BLH(ol, ustar, lat, zm):
    # --- if no boundary layer height available use Kljun et al., 2015 analytical solution for Habl
    #   blh for stable and neutral conditions - Nieuwstadt (1981)
    #            h = (L/3.8)*(-1 + (1 + 2.28*(ustar/(f*L)))^0.5)                          (Eq.1.1)
    #              with L = MOL, ustar - friction velocity, g = acceleration due to gravity and
    #                   f = 2 omega sin (phi) = Coriolis parameter
    #   blh for convective conditions (needs to be integrated due to near symmetric diurnal cycle
    #            of surface heat flux). The resulting rate of change for the boundary layer height
    #            is implicit in h and may be solved iteratively or using h(t_i) to determine the rate
    #            of change to yield h(t_i+1)
    #            at sunrise before Fh becomes positive - use Eq.1.1 for initial conditions
    #               then dh/dt = bar(w'Theta')_0 / gamma*[(h^2 / ((1+2*A)*h - 2*B*k*L))   (Eq.1.2)
    #                           + ((C * ustar^2 * T) / (gamma*g*[(1+A)*h - B*k*L]))]^(-1)
    #            gamma = gradient of potential temp above convective bl; ~0.01 K/m for typical midlatitude
    #                      (dh/dt quite sensitive to gamma)
    #            A = 0.2; B + 2.5; C = 8 (derived from similarity relations in cbl)
    #   so make sure H is availble from external sources, you really don't want to calculate it
    #A = 0.2
    #B = 2.5
    #C = 8
    #gamma = 0.01          # K/m
    omega = 0.000072921   # rad/s
    f = 2.0 * omega * numpy.sin(lat * numpy.pi / 180.)
    if ol > 0:    # stable
        Habl = (ol/3.8)*(-1 + 2.28*(ustar/(f*ol)))**0.5
    elif ol <= 0: # convective
        # here we need to model the data day by day
        Habl = (ol/3.8)*(-1 + 2.28*(ustar/(f*ol)))**0.5
    elif float(zm)/ol <= -15.5:
        # not valid in Kljun et al., 2015
        Habl = -9999
    return Habl

def create_index_list(cf, d, date):
    """
    Create a list of indices of the datetime for the requested climatology
    Single  = only one element for Start and End, here range is index , index+1
            difficulty, the time for this single element is actually the timestep
            from previous to the named in the list (timestamp at end of interval
    Special = only one element for Start and End
    Daily   = get a climatology for each day, forget about the couple of hours
            before and after the first day, may be able to use
            qcutils.GetDateIndex(ldt,date_str,ts=30,default=0,match='startnextday')
            qcutils.GetDateIndex(ldt,date_str,ts=30,default=0,match='endpreviousday')
    Monthly =
    Annual  =
    """
    #import datetime
    climfreq = d["Climatology"]
    # firstly what climatology is requested
    if climfreq == 'Single':
        list_StDate = [None]*1
        list_EnDate = [None]*1
        if 'StartDate' in cf['Options'].keys():
            xlStDate = cf['Options']['StartDate']
            list_StDate[0] = qcutils.GetDateIndex(date,xlStDate,ts=d["flux_period"],default=0,match='exact')
        else:
            logger.error("No StartDate given. Define which time for footprint calculation in StartDate (DD/MM/YYYY hh:mm)")

        list_EnDate[0] = list_StDate[0]+1

    elif climfreq == 'Special':
        list_StDate = []
        list_EnDate = []
        if 'StartDate' in cf['Options'].keys():
            xlStDate = cf['Options']['StartDate']
            list_StDate[0] = qcutils.GetDateIndex(date,xlStDate,ts=d["flux_period"],default=0,match='exact')
        else:
            list_StDate[0]  = 0               # start from begin of file
        if 'EndDate' in cf['Options'].keys():
            xlEnDate = cf['Options']['EndDate']
            list_EnDate[0] = qcutils.GetDateIndex(date,xlEnDate,ts=d["flux_period"],default=0,match='exact')
        else:
            list_EnDate[0]  = len(date)-1 # run to end of file

    elif climfreq == 'Daily':
        StDate = date[0]
        EnDate = date[-1]
        sd = pandas.date_range(start=StDate, end=EnDate, freq='D', normalize=True)    # frequency daily
        ndays     = len(sd)
        list_StDate = []
        list_EnDate = []
        list_StDate.append(qcutils.GetDateIndex(date,sd[0],ts=d["flux_period"],default=0,match='exact'))
        list_EnDate.append(qcutils.GetDateIndex(date,sd[1],ts=d["flux_period"],default=-1,match='exact'))
        for i in range(1,ndays-1):
            list_StDate.append(qcutils.GetDateIndex(date,sd[i],ts=d["flux_period"],default=0,match='exact') +1)
            list_EnDate.append(qcutils.GetDateIndex(date,sd[i+1],ts=d["flux_period"],default=-1,match='exact'))
        test_i = qcutils.GetDateIndex(date,sd[-1],ts=d["flux_period"],default=0,match='exact')
        if test_i < len(date)-2: # at least one value for the next day, so only midnight not allowed
            list_StDate.append(test_i+1)
            list_EnDate.append(len(date)-1)

    elif climfreq == 'Monthly':
        StDate = date[0]
        EnDate = date[-1]
        sm = pandas.date_range(start=StDate, end=EnDate, freq='MS', normalize=True)    # frequency monthly
        num_int = len(sm)
        list_StDate = []
        list_EnDate = []
        test_i = qcutils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact')
        if test_i > 0:
            list_StDate.append(0)
            list_EnDate.append(test_i)
            list_StDate.append(qcutils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact')+1)
            list_EnDate.append(qcutils.GetDateIndex(date,sm[1],ts=d["flux_period"],default=-1,match='exact'))
        else:
            list_StDate.append(qcutils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact'))
            list_EnDate.append(qcutils.GetDateIndex(date,sm[1],ts=d["flux_period"],default=-1,match='exact'))
        for i in range(1,num_int-1):
            list_StDate.append(qcutils.GetDateIndex(date,sm[i],ts=d["flux_period"],default=0,match='exact')+1)
            list_EnDate.append(qcutils.GetDateIndex(date,sm[i+1],ts=d["flux_period"],default=-1,match='exact'))
        test_i = qcutils.GetDateIndex(date,sm[-1],ts=d["flux_period"],default=0,match='exact')
        if test_i < len(date)-2: # at least one value for the next day, so only midnight not allowed
            list_StDate.append(test_i+1)
            list_EnDate.append(len(date)-1)

    elif climfreq == 'Annual':
        # Find number of years in df
        StDate = date[0]
        EnDate = date[-1]
        years_index = []
        year = date.apply(lambda x: x.year)
        for i in range(min(year),max(year)+1):
            years_index.append(i)
        num = len(years_index)
        years_index.append(max(years_index)+1)
        print num,years_index
        list_StDate = []
        list_EnDate = []
        st = datetime.datetime(years_index[0],1,1,0,0)
        en = datetime.datetime(years_index[1],1,1,0,0)
        list_StDate.append(qcutils.GetDateIndex(date,st,ts=d["flux_period"],default=0,match='exact'))
        list_EnDate.append(qcutils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact'))
        if num > 1:
            if num > 2:
                for i in range(1,num-1):
                    st = datetime.datetime(years_index[i],1,1,0,0)
                    en = datetime.datetime(years_index[i+1],1,1,0,0)
                    list_StDate.append(qcutils.GetDateIndex(date,st,ts=d["flux_period"],default=0,match='exact')+1)
                    list_EnDate.append(qcutils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact'))
            st = datetime.datetime(years_index[num-1],1,1,0,0)
            en = datetime.datetime(years_index[num],1,1,0,0)
            test_is = qcutils.GetDateIndex(date,st,ts=d["flux_period"],default=-1,match='exact')
            test_ie = qcutils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact')
            if test_ie - test_is > 2:
                list_StDate.append(test_is+1)
                list_EnDate.append(test_ie)

    print list_StDate,list_EnDate

    return list_StDate,list_EnDate

#def fp_data_in_original(cf,mode):

    ## read input data and prepare for input into Kormann and Meixner, 2001 or Kljun et al., 2015
    ## python routines
    ## ---------------------- Get input / output file name ------------------------------------

    ## Set input file and output path and create directories for plots and results
    #filepath = cf['Files']['file_path']
    #file_in = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'])
    ## get a dictionary of the variable names
    #var_list = cf["InVariables"].keys()
    ##print 'Variables=',var_list
    #names = {}

    #for item in var_list:
        #if "ncname" in cf["InVariables"][item].keys():
            #names[item] = cf["InVariables"][item]["ncname"]
        #else:
            #names[item] = item
    #ext_var_list = cf["ExtVariables"].keys()
    #for item in ext_var_list:
        #if "ncname" in cf["ExtVariables"][item].keys():
            #names[item] = cf["ExtVariables"][item]["ncname"]
        #else:
            #names[item] = item
    ## add the xlDateTime
    #names["xlDateTime"] = "xlDateTime"
    #names["DateTime"] = "DateTime"
    #names["Day"] = "Day"
    #names["Month"] = "Month"
    #names["Year"] = "Year"
    ## read the netcdf file
    #logger.info(' Reading netCDF file '+file_in)
    #ds = qcio.nc_read_series(file_in)
    #dates_list = ds.series["DateTime"]["Data"]
    #nrecs = int(ds.globalattributes["nc_nrecs"])
    ## read the external file for Habl if mode = kljun
    ## === check to see if we have Habl timeseries in imports ??? What if not? Botheration!
    #if mode=="kljun":
        #qcgf.ImportSeries(cf,ds)
    ## ds now includes all the data, also Habl
    ##for item in ["L","V_Sd","z0"]:
    #if "L" not in ds.series.keys():
        #Ta  = qcutils.GetVariable(ds,"Ta")
        #Ah = qcutils.GetVariable(ds,"Ah")
        #ps = qcutils.GetVariable(ds,"ps")
        #ustar = qcutils.GetVariable(ds,"ustar")
        #Fh = qcutils.GetVariable(ds,"Fh")
    #qcutils.CreateVariable(ds,"L")

    ## now get the data
    #d = {}
    #f = {}

    #for item in names.keys():
        #data,flag,attr = qcutils.GetSeries(ds,names[item])
        ##print item
        #if item == "StdV":
            #if names[item] == "uyuy":
                #data = numpy.sqrt(uyuy)
        #d[item] = data
        #f[item] = flag

    ## if the cross wind standard deviation is not in the data set (quite common) then use something else
    #if "StdV" not in var_list:
        #names["StdV"] = "StdV"
        #d["StdV"] = 0.5*d["Ws"].astype(float)
        #logger.warning("Stdev of cross wind component not in data structure, will be estimated from 0.5*Ws")

    ## === Monin Obukhov length
    #if "L" not in var_list:
        #names["L"] = "L"
        #sTa = d["Ta"].astype(float)
        #sAh = d["Ah"].astype(float)
        #sps = d["ps"].astype(float)
        #sus = d["ustar"].astype(float)
        #sFh = d["Fh"].astype(float)
        ##data = mf.molen(sTa,sAh,sps,sus,sFh,fluxtype="sensible")
        #data = mf.molen(sTa,sAh,sps,sus,sFh,fluxtype="sensible")
        ##print 'calculate L ', sTa[ngive],sAh[ngive],sps[ngive],sus[ngive],sFh[ngive],data[ngive]
        #d["L"] = numpy.where(data==c.missing_value,'None',data)
        ##d["L"] = numpy.where(data==c.missing_value,numpy.nan,data)
        ##d["L"] = data.where(data.notnull(), None)
        #d["L"] = d["L"].astype(float)
    ## === roughness length
    #if "z0" not in var_list:
        #names["z0"] = "z0"
        #zT = float(cf["Tower"]["tower_height"])
        #zC = float(cf["Tower"]["canopy_height"])
        #zm = zT-(2.0/3.0)*zC
        #sL =  d["L"].astype(float)
        #sWs = d["Ws"].astype(float)
        ##data = 0.15* zC # Value as used in EddyPro if z0 not entered as a variable directly
        ##data = z0calc(zm,sL,sWs,sus)
        #data = z0calc(zm,sL,sWs,sus)
        ##print 'calculate z0 ', zm,sL,sWs,sus,temp,data
        ## z0calc(sL,zm,sWs,sus); estimate for z0 using MO-Similarity Theory
        #d["z0"] = numpy.where(data==c.missing_value,'None',data)
        #d["z0"] = d["z0"].astype(float)
    ## === atmospheric boundary layer
    ##if "Habl" not in ext_var_list:
    ##   # see Kljun et al., 2015 paper for an estimate
    ##   logger.warning("Habl boundary layer height not in data structure, will be estimated")
    ##   names["Habl"] = "Habl"
    ##   v1 =
    ##   v2 =
    ##   v3 =
    ##
    ##   d["Habl"] = ????,

    #df=pandas.DataFrame(d,index=dates_list)
    ## replace missing values with NaN
    ##df.replace(c.missing_value,numpy.nan)
    ## Build dictionary of additional configs
    #d={}
    ## === Which climatology, either definded time, daily, monthly or annual
    #d["Climatology"] = qcutils.get_keyvaluefromcf(cf,["Options"],"Climatology",default="Special")
    #climfreq = d["Climatology"]
    ##
    #if "out_filename" in cf['Files']:
        #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['out_filename'])
    #else:
        #if climfreq == 'Annual':
            #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_y_fp.nc"))
        #elif climfreq == 'Monthly':
            #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_m_fp.nc"))
        #elif climfreq == 'Daily':
            #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_d_fp.nc"))
        #elif climfreq == 'Single':
            #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_s_fp.nc"))
        #elif climfreq == 'Special':
            #file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_s_fp.nc"))

    #plot_path = "plots/"
    #if "plot_path" in cf["Files"]: plot_path = os.path.join(cf["Files"]["plot_path"],"FP/")
    #if not os.path.isdir(plot_path): os.makedirs(plot_path)

    #results_path = filepath
    #if not os.path.isdir(results_path): os.makedirs(results_path)

    #d["tower_height"]   = float(cf["Tower"]["tower_height"])
    #d["canopy_height"]  = float(cf["Tower"]["canopy_height"])
    #d["footprint_size"] = int(cf["Tower"]["footprint_size"])
    #d["zm_d"]           = d["tower_height"]-(2.0/3.0*d["canopy_height"])
    #d["xTower"]         = 0 #int(cf['Tower']['xTower'])
    #d["yTower"]         = 0 #int(cf['Tower']['yTower'])
    #d["xmin"]           = -0.5*d["footprint_size"]
    #d["xmax"]           =  0.5*d["footprint_size"]
    #d["ymin"]           = -0.5*d["footprint_size"]
    #d["ymax"]           =  0.5*d["footprint_size"]
    #d["nx"]             = int(cf["Tower"]["num_cells"])

    #d["flux_period"] = int(ds.globalattributes["time_step"])
    ##d["timezone"] = int(ds.globalattributes["timezone"])
    #d["site_name"]   = ds.globalattributes["site_name"]
    #if "Latitude" in cf["Tower"]:
        #d["latitude"] = cf["Tower"]["Latitude"]
    #else:
        #d["latitude"] = ds.globalattributes["latitude"]
    #if "Longitude" in cf["Tower"]:
        #d["longitude"] = cf["Tower"]["Longitude"]
    #else:
        #d["longitude"]   = ds.globalattributes["longitude"]

    #d["call_mode"]   = qcutils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive",mode="quiet")
    #d["show_plots"]  = qcutils.get_keyvaluefromcf(cf,["Options"],"show_plots",default=True,mode="quiet")
    #d["file_out"]    = file_out
    #d["plot_path"]    = plot_path

    #return df, d, ds

def get_footprint_data_in(cf, mode):
    # read input data and prepare for input into Kormann and Meixner, 2001 or Kljun et al., 2015
    # python routines
    # ---------------------- Get input / output file name ------------------------------------

    # Set input file and output path and create directories for plots and results
    file_in = os.path.join(cf['Files']['file_path'], cf['Files']['in_filename'])
    # read the netcdf file
    msg = ' Reading netCDF file ' + str(file_in)
    logger.info(msg)
    ds = qcio.nc_read_series(file_in)
    nrecs = int(ds.globalattributes["nc_nrecs"])
    # array of 0s for QC flag
    f0 = numpy.zeros(nrecs, dtype=numpy.int32)
    # array of 1s for QC flag
    f1 = numpy.ones(nrecs, dtype=numpy.int32)
    # read the external file for Habl if mode = kljun
    # === check to see if we have Habl timeseries in imports ??? What if not? Botheration!
    if mode == "kljun":
        qcgf.ImportSeries(cf, ds)
    else:
        Habl = qcutils.create_empty_variable("Habl", nrecs)
        Habl["Label"] = "Habl"
        Habl["Data"] = numpy.ma.array(numpy.full(nrecs, 1000))
        Habl["Flag"] = f0
        Habl["Attr"] = {"long_name":" Boundary-layer height", "units":"m",
                        "standard_name":"not defined"}
        qcutils.CreateVariable(ds, Habl)
    # cjeck to see if Monin-Obukhov length is in the data structure
    if "L" not in ds.series.keys():
        # if not, calculate it
        qcts.CalculateMoninObukhovLength(ds)
    # if the cross wind standard deviation is not in the data set (quite common) then use something else
    if "V_Sd" not in ds.series.keys():
        # could do better with:
        # 1) reprocess L3 and output variance of U, V and W
        # 2) estimated from standard deviation of wind direction (if available)
        # 3) estimate using MO relations (needs Habl)
        logger.warning("Stdev of cross wind component not in data structure, estimated as 0.5*Ws")
        Ws = qcutils.GetVariable(ds, "Ws")
        V_Sd = qcutils.create_empty_variable("V_Sd", nrecs)
        V_Sd["Data"] = 0.5*Ws["Data"]
        V_Sd["Flag"] = numpy.where(numpy.ma.getmaskarray(V_Sd["Data"])==True, f1, f0)
        V_Sd["Attr"]["long_name"] = "Variance of cross-wind velocity component, estimated from Ws"
        V_Sd["Attr"]["units"] = "(m/s)2"
        V_Sd["Attr"]["height"] = Ws["Attr"]["height"]
        qcutils.CreateVariable(ds, V_Sd)
    # === roughness length
    if "z0" not in ds.series.keys():
        z0 = qcutils.create_empty_variable("z0", nrecs)
        # check the global attriibutes first
        if "roughness_length" in ds.globalattributes.keys():
            roughness_length = float(ds.globalattributes["roughness_length"])
            z0["Data"] = numpy.ma.array(numpy.full(nrecs, roughness_length))
            z0["Attr"]["long_name"] = "Roughness length from global attributes"
        elif "roughness_length" in cf["Tower"]:
            roughness_length = float(cf["Tower"]["roughness_length"])
            z0["Data"] = numpy.ma.array(numpy.full(nrecs, roughness_length))
            z0["Attr"]["long_name"] = "Roughness length from footprint control file"
        else:
            zT = float(cf["Tower"]["tower_height"])
            zC = float(cf["Tower"]["canopy_height"])
            zm = zT-(2.0/3.0)*zC
            L = qcutils.GetVariable(ds, "L")
            ustar = qcutils.GetVariable(ds, "ustar")
            Ws = qcutils.GetVariable(ds, "Ws")
            z0["Data"] = z0calc(zm, L["Data"], Ws["Data"], ustar["Data"])
            z0["Attr"]["long_name"] = "Roughness length calculated from u*, L, Ws and (z-d)"
        z0["Flag"] = numpy.where(numpy.ma.getmaskarray(z0["Data"])==True, f1, f0)
        z0["Attr"]["units"] = "m"
        qcutils.CreateVariable(ds, z0)
    return ds

def get_footprint_cfg(cf, ds):
    # Build dictionary of additional configs
    d={}
    # === Which climatology, either definded time, daily, monthly or annual
    d["Climatology"] = qcutils.get_keyvaluefromcf(cf,["Options"],"Climatology",default="Special")
    climfreq = d["Climatology"]
    #
    if "out_filename" in cf['Files']:
        file_out = os.path.join(cf['Files']['file_path'],cf['Files']['out_filename'])
    else:
        if climfreq == 'Annual':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_y_fp.nc"))
        elif climfreq == 'Monthly':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_m_fp.nc"))
        elif climfreq == 'Daily':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_d_fp.nc"))
        elif climfreq == 'Single':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_s_fp.nc"))
        elif climfreq == 'Special':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_s_fp.nc"))

    plot_path = "plots/"
    if "plot_path" in cf["Files"]: plot_path = os.path.join(cf["Files"]["plot_path"],"FP/")
    if not os.path.isdir(plot_path): os.makedirs(plot_path)

    results_path = cf['Files']['file_path']
    if not os.path.isdir(results_path): os.makedirs(results_path)

    d["tower_height"]   = float(cf["Tower"]["tower_height"])
    d["canopy_height"]  = float(cf["Tower"]["canopy_height"])
    d["footprint_size"] = int(cf["Tower"]["footprint_size"])
    d["zm_d"]           = d["tower_height"]-(2.0/3.0*d["canopy_height"])
    d["xTower"]         = 0 #int(cf['Tower']['xTower'])
    d["yTower"]         = 0 #int(cf['Tower']['yTower'])
    d["xmin"]           = -0.5*d["footprint_size"]
    d["xmax"]           =  0.5*d["footprint_size"]
    d["ymin"]           = -0.5*d["footprint_size"]
    d["ymax"]           =  0.5*d["footprint_size"]
    d["nx"]             = int(cf["Tower"]["num_cells"])

    d["flux_period"] = int(ds.globalattributes["time_step"])
    #d["timezone"] = int(ds.globalattributes["timezone"])
    d["site_name"] = ds.globalattributes["site_name"]
    if "Latitude" in cf["Tower"]:
        d["latitude"] = cf["Tower"]["Latitude"]
    else:
        d["latitude"] = ds.globalattributes["latitude"]
    if "Longitude" in cf["Tower"]:
        d["longitude"] = cf["Tower"]["Longitude"]
    else:
        d["longitude"] = ds.globalattributes["longitude"]

    d["call_mode"] = qcutils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive",mode="quiet")
    d["show_plots"] = qcutils.get_keyvaluefromcf(cf,["Options"],"show_plots",default=True,mode="quiet")
    d["file_out"] = file_out
    d["plot_path"] = plot_path

    return d