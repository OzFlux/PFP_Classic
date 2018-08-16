# standard modules
import csv
import datetime
import logging
import os
import sys
# 3rd party modules
#import dateutil
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile
import netCDF4
import numpy
import pandas
import pysolar
import xlrd
# PFP modules
import constants as c
import meteorologicalfunctions as mf
import pfp_io
import pfp_utils

from windrose import WindroseAxes


logger = logging.getLogger("pfp_log")

if not os.path.exists("./scripts/"):
    print "PyFluxPro: the scripts directory is missing"
    sys.exit()
sys.path.append('scripts')

c_d2r=c.Pi/180.0

# Coordinate steps in footprint process
def windrose_main(cf):
    """
    This script reads data from a PyFluxPro .nc file and processes the data for a
         windrose climatology
    === > output for the climatology
    (a) daily climatology
    (b) monthly climatology
    (c) seasonal climatology
    (d) annual climatology
    GOAL: windrose plot for either annual, seasonal monthly or special
    DONE: copied footprint script and adjusted it to deal with windrose data, 
          kept annual, monthly, special time setting
          copied windrose.py into script directory
    
    C.M.Ewenz, 12 June 2018
    Latest changes: 6 July 2018 - fixed the climatology time axis
    """
    logger.info(' Read input data files ...')
    # get the L3 data
    ds = get_windrose_data_in(cf)
    ldt = ds.series["DateTime"]["Data"]
    # get the configuration data for the footprint
    d = get_windrose_cfg(cf, ds)
    logger.info(' Starting windrose climatology ...')
    list_StDate, list_EnDate = create_index_list(cf, d, ldt)
    # After deciding which climatology is done, let's do it!
    for i in range(0,len(list_StDate)):
        # get the start and end indices
        si = list_StDate[i]
        ei = list_EnDate[i]
        # get the series as masked arrays
        wind_speed, _, _ = pfp_utils.GetSeriesasMA(ds, "Ws", si=si, ei=ei)
        wind_dir, _, _   = pfp_utils.GetSeriesasMA(ds, "Wd", si=si, ei=ei)
        # get a composite mask over all variables
        mask_all = numpy.ma.getmaskarray(wind_speed)
        for item in [wind_dir]:
            mask_item = numpy.ma.getmaskarray(item)
            mask_all = numpy.ma.mask_or(mask_all, mask_item)
        # and then apply the composite mask to all variables and remove masked elements
        wind_speed = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, wind_speed)))
        wind_dir = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, wind_dir)))
        if len(wind_speed) == 0:
            msg = "No windrose input data for "+str(ldt[si])+" to "+str(ldt[ei])
            logger.warning(msg)
            num[irun]=0
        else:
            # get the default plot width and height
            PlotWidth  = float(cf['General']['PlotWidth'])
            PlotHeight = float(cf['General']['PlotHeight'])
            # get the image filename if plotting if OzFlux_area_image in cf['General'].keys():
            plotname = 'plots/Windrose ' + d["site_name"].replace(' ','') + str(list_StDate[i]).zfill(6) +'.jpg'
            #Figure 1: windrose like a stacked histogram with normed (displayed in percent) results
            fig = plt.figure(i,figsize=(PlotWidth, PlotHeight), dpi=80, facecolor='w', edgecolor='w')
            rect = [0.1, 0.1, 0.8, 0.8]
            ax = WindroseAxes(fig, rect, facecolor='w')
            fig.add_axes(ax)
            ax.bar(wind_dir, wind_speed, normed=True, opening=0.8, edgecolor='white')
            #set_legend(ax)
            l = ax.legend()
            #l = ax.legend(axespad=-0.10)
            plt.setp(l.get_texts(), fontsize=8)
            
            ntime1 = ldt[si].strftime("%Y%m%d%H%M")
            ntime2 = ldt[ei].strftime("%Y%m%d%H%M")
            if d["Climatology"] == 'Single':
                text = 'Windrose ' + d["site_name"] + ' ' + ntime1
            else:
                text = 'Windrose ' + d["site_name"] + ' ' + ntime1 + '  to  ' + ntime2
            plt.title(text)
            
            plt.savefig(plotname)
            
            fig.show()
        
    logger.info(" Finished windrose plotting")
 
    # ================================================================
def create_index_list(cf, d, date):
    """
    Create a list of indices of the datetime for the requested climatology
    Single  = only one element for Start and End, here range is index , index+1
            difficulty, the time for this single element is actually the timestep
            from previous to the named in the list (timestamp at end of interval
    Special = only one element for Start and End
    Daily   = get a climatology for each day, forget about the couple of hours
            before and after the first day, use
            pfp_utils.GetDateIndex(ldt,date_str,ts=30,default=0,match='startnextday')
            pfp_utils.GetDateIndex(ldt,date_str,ts=30,default=0,match='endpreviousday')
    Monthly =
    Annual  =
    """
    #import datetime
    climfreq = d["Climatology"]
    # firstly what climatology is requested
    if climfreq == 'Single':
        list_StDate = []
        list_EnDate = []
        if 'StartDate' in cf['Options'].keys():
            xlStDate = cf['Options']['StartDate']
            list_StDate.append(pfp_utils.GetDateIndex(date,xlStDate,ts=d["flux_period"],default=0,match='exact'))
        else:
            logger.error("No StartDate given. Define which time for footprint calculation in StartDate (DD/MM/YYYY hh:mm)")

        list_EnDate.append(list_StDate[0]+1)

    elif climfreq == 'Special':
        list_StDate = []
        list_EnDate = []
        if 'StartDate' in cf['Options'].keys():
            xlStDate = cf['Options']['StartDate']
            print xlStDate
            list_StDate.append(pfp_utils.GetDateIndex(date,xlStDate,ts=d["flux_period"],default=0,match='exact'))
        else:
            list_StDate.append(0)         # start from begin of file
        if 'EndDate' in cf['Options'].keys():
            xlEnDate = cf['Options']['EndDate']
            list_EnDate.append(pfp_utils.GetDateIndex(date,xlEnDate,ts=d["flux_period"],default=0,match='exact'))
        else:
            list_EnDate.append(len(date)-1) # run to end of file

    elif climfreq == 'Daily':
        StDate = date[0]
        EnDate = date[-1]
        sd = pandas.date_range(start=StDate, end=EnDate, freq='D', normalize=True)    # frequency daily
        ndays     = len(sd)
        list_StDate = []
        list_EnDate = []
        list_StDate.append(pfp_utils.GetDateIndex(date,sd[0],ts=d["flux_period"],default=0,match='exact'))
        list_EnDate.append(pfp_utils.GetDateIndex(date,sd[1],ts=d["flux_period"],default=-1,match='exact'))
        for i in range(1,ndays-1):
            list_StDate.append(pfp_utils.GetDateIndex(date,sd[i],ts=d["flux_period"],default=0,match='exact') +1)
            list_EnDate.append(pfp_utils.GetDateIndex(date,sd[i+1],ts=d["flux_period"],default=-1,match='exact'))
        test_i = pfp_utils.GetDateIndex(date,sd[-1],ts=d["flux_period"],default=0,match='exact')
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
        test_i = pfp_utils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact')
        if test_i > 0:
            list_StDate.append(0)
            list_EnDate.append(test_i)
            list_StDate.append(pfp_utils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact')+1)
            list_EnDate.append(pfp_utils.GetDateIndex(date,sm[1],ts=d["flux_period"],default=-1,match='exact'))
        else:
            list_StDate.append(pfp_utils.GetDateIndex(date,sm[0],ts=d["flux_period"],default=0,match='exact'))
            list_EnDate.append(pfp_utils.GetDateIndex(date,sm[1],ts=d["flux_period"],default=-1,match='exact'))
        for i in range(1,num_int-1):
            list_StDate.append(pfp_utils.GetDateIndex(date,sm[i],ts=d["flux_period"],default=0,match='exact')+1)
            list_EnDate.append(pfp_utils.GetDateIndex(date,sm[i+1],ts=d["flux_period"],default=-1,match='exact'))
        test_i = pfp_utils.GetDateIndex(date,sm[-1],ts=d["flux_period"],default=0,match='exact')
        if test_i < len(date)-2: # at least one value for the next day, so only midnight not allowed
            list_StDate.append(test_i+1)
            list_EnDate.append(len(date)-1)

    elif climfreq == 'Annual':
        # Find number of years in df
        StDate = date[0]
        EnDate = date[-1]
        years_index = []
        #date.apply(lambda x: x.year)
        #for i in range(min(year),max(year)+1):
        for i in range(StDate.year,EnDate.year+1):
            years_index.append(i)
        num = len(years_index)
        years_index.append(max(years_index)+1)
        #print num,years_index
        list_StDate = []
        list_EnDate = []
        st = datetime.datetime(years_index[0],1,1,0,0)
        en = datetime.datetime(years_index[1],1,1,0,0)
        list_StDate.append(pfp_utils.GetDateIndex(date,st,ts=d["flux_period"],default=0,match='exact'))
        list_EnDate.append(pfp_utils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact'))
        if num > 1:
            if num > 2:
                for i in range(1,num-1):
                    st = datetime.datetime(years_index[i],1,1,0,0)
                    en = datetime.datetime(years_index[i+1],1,1,0,0)
                    list_StDate.append(pfp_utils.GetDateIndex(date,st,ts=d["flux_period"],default=0,match='exact')+1)
                    list_EnDate.append(pfp_utils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact'))
            st = datetime.datetime(years_index[num-1],1,1,0,0)
            en = datetime.datetime(years_index[num],1,1,0,0)
            test_is = pfp_utils.GetDateIndex(date,st,ts=d["flux_period"],default=-1,match='exact')
            test_ie = pfp_utils.GetDateIndex(date,en,ts=d["flux_period"],default=-1,match='exact')
            if test_ie - test_is > 2:
                list_StDate.append(test_is+1)
                list_EnDate.append(test_ie)
    return list_StDate,list_EnDate

def get_windrose_data_in(cf):

    # read input data and prepare for input into Kljun et al., 2015 python routines
    # ---------------------- Get input / output file name ------------------------------------
    
    # Set input file and output path and create directories for plots and results
    filepath = cf['Files']['file_path']
    #file_in = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'])
    #
    # read input data and prepare for input into Kormann and Meixner, 2001 or Kljun et al., 2015
    # python routines
    # ---------------------- Get input / output file name ------------------------------------

    # Set input file and output path and create directories for plots and results
    file_in = os.path.join(cf['Files']['file_path'], cf['Files']['in_filename'])
    # read the netcdf file
    msg = ' Reading netCDF file ' + str(file_in)
    logger.info(msg)
    ds = pfp_io.nc_read_series(file_in)
    nrecs = int(ds.globalattributes["nc_nrecs"])
    # array of 0s for QC flag
    f0 = numpy.zeros(nrecs, dtype=numpy.int32)
    # array of 1s for QC flag
    f1 = numpy.ones(nrecs, dtype=numpy.int32)

    return ds

def get_windrose_cfg(cf,ds):
    # Build dictionary of additional configs
    plot_path = "plots/"
    if "plot_path" in cf["Files"]: plot_path = os.path.join(cf["Files"]["plot_path"],"WR/")
    if not os.path.isdir(plot_path): os.makedirs(plot_path)

    #results_path = filepath
    #if not os.path.isdir(results_path): os.makedirs(results_path)

    d={}
    # === Which climatology, either definded time, daily, monthly or annual
    d["Climatology"] = pfp_utils.get_keyvaluefromcf(cf,["Options"],"Climatology",default="Special")
    climfreq = d["Climatology"]

    d["tower_height"]   =float(cf["Tower"]["tower_height"])
    d["canopy_height"]  =float(cf["Tower"]["canopy_height"])
    d["footprint_size"] =int(cf["Tower"]["footprint_size"])
    d["zm_d"]           = d["tower_height"]-(2.0/3.0*d["canopy_height"])
    d["xTower"]         = 0 #int(cf['Tower']['xTower'])
    d["yTower"]         = 0 #int(cf['Tower']['yTower'])
    d["xmin"]           = -0.5*d["footprint_size"]
    d["xmax"]           =  0.5*d["footprint_size"]
    d["ymin"]           = -0.5*d["footprint_size"]
    d["ymax"]           =  0.5*d["footprint_size"]
    d["nx"]             = int(cf["Tower"]["num_cells"])
    #d["start_time"]     = dates_list[0]
    #d["end_time"]       = dates_list[-1]
    

    d["flux_period"] = int(ds.globalattributes["time_step"])
    d["site_name"]   = ds.globalattributes["site_name"]

    if "Latitude" in cf["Tower"]:
        d["latitude"] = cf["Tower"]["Latitude"]
    else:
        d["latitude"] = ds.globalattributes["latitude"]
    if "Longitude" in cf["Tower"]:
        d["longitude"] = cf["Tower"]["Longitude"]
    else:
        d["longitude"] = ds.globalattributes["longitude"]

    d["call_mode"]   = pfp_utils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive",mode="quiet")
    d["show_plots"]  = pfp_utils.get_keyvaluefromcf(cf,["Options"],"show_plots",default=True,mode="quiet")
    d["plot_out"]    = plot_path

    return d

def xldate_to_datetime(xldate):
    import datetime
    tempDate = datetime.datetime(1899,12,30)
    deltaDays = datetime.timedelta(days=int(xldate))
    secs = (int((xldate%1)*86400))
    deltaSeconds = datetime.timedelta(seconds=secs)
    TheTime = (tempDate + deltaDays + deltaSeconds )
    return TheTime.strftime("%d-%m-%Y %H:%M")
