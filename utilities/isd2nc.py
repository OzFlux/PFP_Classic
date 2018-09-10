#%run basics
# Python modules
import calendar
from collections import OrderedDict
import copy
import csv
import datetime
import gzip
import logging
import os
import pytz
import sys
import time
# 3rd party modules
from configobj import ConfigObj
import matplotlib.pyplot as plt
import netCDF4
import numpy
import scipy
import scipy.stats
from scipy.interpolate import InterpolatedUnivariateSpline
import xlrd
import xlwt
# pfp modules
if not os.path.exists("../scripts/"):
    print "PyFluxPro: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('../scripts')
import constants as c
import meteorologicalfunctions as mf
import pfp_io
import pfp_log
import pfp_utils

t = time.localtime()
rundatetime = datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]).strftime("%Y%m%d%H%M")
log_filename = 'isd2nc_'+rundatetime+'.log'
logger = pfp_log.init_logger(logger_name="pfp_log", file_handler=log_filename)

def average_duplicate_times(ds_in, time_step):
    """
    Purpose:
     Remove duplicate time steps by averaging data with the same time stamp.
     The routine uses scipy.stats.binned_statistics() to bin the data based
     on the time (bins have width time_step and are centered on times that
     are an integral of time_step).
    Usage:
     ds_out = average_duplicate_times(ds_in, time_step=30)
    Side effects:
     The time given for the averages and sums is the end of the time period.
    Author: PRI
    Date: October 2017
    """
    logger.info("Getting data onto a regular time step")
    # get the time as a number (see attr["units"] for units)
    time_var = pfp_utils.GetVariable(ds_in, "time")
    # generate an array of bin edges for use by binned_statistics()
    bin_width = time_step*60
    # round the ISD start time to an integral of the time step
    t0 = time_step*60*int(time_var["Data"].data[0]/(time_step*60))
    bin_first = t0 - bin_width
    # round the ISD end time to an integral of the time step
    t1 = time_step*60*int(time_var["Data"].data[-1]/(time_step*60))
    # make sure we go 1 beyond the end time
    if t1 < time_var["Data"][-1]:
        t1 = t1 + bin_width
    # generate an array of bin edges
    bin_last = t1 + bin_width
    bins = numpy.arange(bin_first, bin_last, bin_width)
    # get the number of records in the output series
    nrecs = len(bins)-1
    # generate series of zeros and ones to be used as QC flags
    f0 = numpy.zeros(nrecs)
    f1 = numpy.ones(nrecs)
    # create an output data structure with a copy of the input global attributes
    ds_out = pfp_io.DataStructure()
    ds_out.globalattributes = copy.deepcopy(ds_in.globalattributes)
    # update the number of records
    ds_out.globalattributes["nc_nrecs"] = nrecs
    # get a list of variable labels but exclude the datetime, time, wind speed and direction variables
    # NB: Wind velocity components U and V will be averaged and wind speed and direction calculated
    # from these.
    labels = [label for label in ds_in.series.keys() if label not in ["DateTime", "time"]]
    # loop over variables
    for label in labels:
        # get the variable
        var_in = pfp_utils.GetVariable(ds_in, label)
        # indices of non-masked elements
        idx = numpy.ma.where(numpy.ma.getmaskarray(var_in["Data"]) == False)[0]
        # check to see if we have at least 1 data point to deal with
        if len(idx) != 0:
            # get the non-masked data as an ndarray
            data_in = numpy.array(var_in["Data"][idx].data)
            time_in = numpy.array(time_var["Data"][idx].data)
            # use binned_statistic() to average records with the same datetime
            if var_in["Label"][0:1] == "P" and var_in["Attr"]["units"] in ["m", "mm"]:
                # do sum for precipitation
                sums, edges, indices = scipy.stats.binned_statistic(time_in, data_in, statistic="sum", bins=bins)
                # convert output to a masked array and mask empty bins
                data_out = numpy.ma.masked_where(numpy.isfinite(sums) == False, numpy.ma.array(sums))
            else:
                # do average for everything else
                means, edges, indices = scipy.stats.binned_statistic(time_in, data_in, statistic="mean", bins=bins)
                # convert output to a masked array and mask empty bins
                data_out = numpy.ma.masked_where(numpy.isfinite(means) == False, numpy.ma.array(means))
            # generate the QC flag
            flag_out = numpy.where(numpy.ma.getmaskarray(data_out) == True, f1, f0)
            # and create the output variable
            var_out = {"Label":label, "Data":data_out, "Flag":flag_out, "Attr":var_in["Attr"]}
        else:
            # no data, so create an empty output variable
            var_out = {"Label":label, "Data":numpy.ma.masked_all(nrecs), "Flag":f1,
                       "Attr":var_in["Attr"]}
        # and write the output variable to the output data structure
        pfp_utils.CreateVariable(ds_out, var_out)
    # generate a series of the bin mid-points
    mids = edges[1:]
    # and convert these to a series of Python datetimes
    attr = copy.deepcopy(ds_in.series["DateTime"]["Attr"])
    ldt_out = {"Label":"DateTime", "Data":netCDF4.num2date(mids, time_var["Attr"]["units"]),
               "Flag":f0, "Attr":attr}
    # and write the datetime to the output data structure
    pfp_utils.CreateVariable(ds_out, ldt_out)
    pfp_utils.get_nctime_from_datetime(ds_out)
    # get wind speed and direction from components
    U = pfp_utils.GetVariable(ds_out, "u")
    V = pfp_utils.GetVariable(ds_out, "v")
    WS, WD = pfp_utils.convert_UVtoWSWD(U, V)
    pfp_utils.CreateVariable(ds_out, WS)
    pfp_utils.CreateVariable(ds_out, WD)
    return ds_out

def contiguous_regions(condition, max_interp_length = 12):
    """
    Purpose:
     Finds contiguous True regions of the boolean array "condition". Returns
     a 2D array where the first column is the start index of the region and the
     second column is the end index.
    Author: Joe Kington (via StackOverflow)
    Date: September 2014
    """
    # Find the indicies of changes in "condition"
    d = numpy.diff(condition)
    idx, = d.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = numpy.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = numpy.r_[idx, condition.size] # Edit
    # Reshape the result into two columns
    idx.shape = (-1,2)
    # filter out contiguous regions whose duration is less than the maximum
    gap_length = idx[:,1] - idx[:,0]
    gl_idx = numpy.where(gap_length > max_interp_length)[0]
    return idx[gl_idx,:]

def convert_time_zone(ds, from_time_zone_in, to_time_zone_in):
    """
    Purpose:
     Convert the datetime series in a data structure from one timezone to another.
    Usage:
    Author: PRI
    Date: June 2017
    """
    from_time_zone = pytz.timezone(from_time_zone_in)
    to_time_zone = pytz.timezone(to_time_zone_in)
    # local pointer to datetime
    from_dt_naive = ds.series["DateTime"]["Data"]
    # local datetime is currently time zone naive, make it aware
    from_dt_aware = [dt.replace(tzinfo=from_time_zone) for dt in from_dt_naive]
    # now the datetime is time zone aware we can convert to another time zone
    to_dt_aware = [dt.astimezone(to_time_zone) for dt in from_dt_aware]
    # take out any day light saving time
    to_dt_aware_nodst = [dt - dt.dst() for dt in to_dt_aware]
    # and make the datetime naive again
    to_dt_naive_nodst = [dt.replace(tzinfo=None) for dt in to_dt_aware_nodst]
    # and write it back to the data structure (shouldn't use direct writes ...)
    ds.series["DateTime"]["Data"] = to_dt_naive_nodst
    ds.series["DateTime"]["Attr"]["time_zone"] = to_time_zone_in
    # and update the time variable
    pfp_utils.get_nctime_from_datetime(ds)
    return

def interpolate_1d(x1, y1, x2):
    """
    Purpose:
     Interpolate data from one time step to another.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    # off we go
    if numpy.ma.is_masked(y1):
        # check we have at least 2 non-masked points
        if numpy.ma.count(y1) >= 2:
            # input Y array is a masked array
            idx = numpy.where(numpy.ma.getmaskarray(y1) == False)[0]
            int_fn = scipy.interpolate.Akima1DInterpolator(x1[idx], y1[idx].data)
            y2 = int_fn(x2)
            y2 = numpy.ma.fix_invalid(y2)
        else:
            msg = "Not enough points (<2) to interpolate"
            logger.warning(msg)
            y2 = numpy.ma.ones(len(x2))*float(c.missing_value)
            y2 = numpy.ma.masked_all_like(y2)
    else:
        int_fn = scipy.interpolate.Akima1DInterpolator(x1, y1)
        y2 = int_fn(x2)
    return y2

def interpolate_ds(ds_in, ts):
    """
    Purpose:
     Interpolate the contents of a data structure onto a different time step.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    logger.info("Interpolating data")
    # instance the output data structure
    ds_out = pfp_io.DataStructure()
    # copy the global attributes
    ds_out.globalattributes = copy.deepcopy(ds_in.globalattributes)
    # add the time step
    ds_out.globalattributes["time_step"] = str(ts)
    # generate a regular time series at the required time step
    dt = ds_in.series["DateTime"]["Data"]
    dt0 = pfp_utils.rounddttots(dt[0], ts=ts)
    if dt0 < dt[0]:
        dt0 = dt0 + datetime.timedelta(minutes=ts)
    dt1 = pfp_utils.rounddttots(dt[-1], ts=ts)
    if dt1 > dt[-1]:
        dt1 = dt1 - datetime.timedelta(minutes=ts)
    idt = [result for result in pfp_utils.perdelta(dt0, dt1, datetime.timedelta(minutes=ts))]
    x1 = numpy.array([toTimestamp(dt[i]) for i in range(len(dt))])
    x2 = numpy.array([toTimestamp(idt[i]) for i in range(len(idt))])
    # loop over the series in the data structure and interpolate
    flag = numpy.zeros(len(idt), dtype=numpy.int32)
    attr = {"long_name":"Datetime", "units":"none"}
    ldt_var = {"Label":"DateTime", "Data":idt, "Flag":flag, "Attr":attr}
    pfp_utils.CreateVariable(ds_out, ldt_var)
    pfp_utils.get_nctime_from_datetime(ds_out)
    nrecs = len(idt)
    ds_out.globalattributes["nc_nrecs"] = nrecs
    # first, we do the air temperature, dew point temperature and surface pressure
    f0 = numpy.zeros(nrecs, dtype=numpy.int32)
    f1 = numpy.ones(nrecs, dtype=numpy.int32)
    for label in ["Ta", "Td", "ps", "RH", "Ah", "q"]:
        var_out = pfp_utils.create_empty_variable(label, nrecs, datetime=idt)
        var_in = pfp_utils.GetVariable(ds_in, label)
        var_out["Data"] = interpolate_1d(x1, var_in["Data"], x2)
        var_out["Flag"] = numpy.where(numpy.ma.getmaskarray(var_out["Data"])==True, f1, f0)
        var_out["Attr"] = copy.deepcopy(var_in["Attr"])
        pfp_utils.CreateVariable(ds_out, var_out)
    # now clamp the dew point so that TD <= TA
    Ta = pfp_utils.GetVariable(ds_out, "Ta")
    Td = pfp_utils.GetVariable(ds_out, "Td")
    Td["Data"] = numpy.ma.where(Td["Data"]<=Ta["Data"], x=Td["Data"], y=Ta["Data"])
    pfp_utils.CreateVariable(ds_out, Td)
    # now we do wind speed and direction by converting to U and V components
    interpolate_wswd(ds_in, x1, ds_out, x2)
    # and lastly, do precipitation
    interpolate_precip(ds_in, x1, ds_out, x2)

    return ds_out

def interpolate_precip(ds_in, x1, ds_out, x2):
    """
    Purpose:
     Transfer precipitation data from one time stamp to another.
     Interpolating precipitation does not make sense since interpolation
     will change the total amount of the precipitation and because
     precipitation is very intermittent, we do not know how it will
     behave at the missing times.  Instead of interpolation, we assign
     all of the precipitation total at the input time stamp to the
     matching output time stamp.
    Assumptions:
     That the input time stamps are a subset of the output time stamps.
    Usage:
    Side effects:
    Author: PRI
    Date: December 2017
    """
    nrecs = ds_out.globalattributes["nc_nrecs"]
    # local pointer to the output datetime
    ldt = pfp_utils.GetVariable(ds_out, "DateTime")
    # just assign any precipitation to the ISD time stamp.
    precip_out = pfp_utils.create_empty_variable("Precip", nrecs, datetime=ldt["Data"])
    # zero the data array
    precip_out["Data"] = precip_out["Data"] * float(0)
    precip_in = pfp_utils.GetVariable(ds_in, "Precip")
    # get indices of elements where the times match
    idx = numpy.searchsorted(x2, numpy.intersect1d(x2, x1))
    precip_out["Data"][idx] = precip_in["Data"]
    precip_out["Flag"] = numpy.zeros(nrecs, dtype=numpy.int32)
    precip_out["Attr"] = copy.deepcopy(precip_in["Attr"])
    pfp_utils.CreateVariable(ds_out, precip_out)
    return

def interpolate_wswd(ds_in, x1, ds_out, x2):
    """
    Purpose:
     Interpolate wind speed and direction by converting them to U and V
     components first, doing the interpolation and then converting back
     to wind speed and direction.
    Usage:
    Side effects
    Author: PRI
    Date: December 2017
    """
    # get the number of records in the output data
    nrecs = ds_out.globalattributes["nc_nrecs"]
    f0 = numpy.zeros(nrecs, dtype=numpy.int32)
    f1 = numpy.ones(nrecs, dtype=numpy.int32)
    # local pointer to the output datetime
    ldt = pfp_utils.GetVariable(ds_out, "DateTime")
    # create empty variables for the output dara
    Ws_out = pfp_utils.create_empty_variable("Ws", nrecs, datetime=ldt["Data"])
    Wd_out = pfp_utils.create_empty_variable("Wd", nrecs, datetime=ldt["Data"])
    U_out = pfp_utils.create_empty_variable("u", nrecs, datetime=ldt["Data"])
    V_out = pfp_utils.create_empty_variable("v", nrecs, datetime=ldt["Data"])
    # get the input wind speed and direction
    Ws_in = pfp_utils.GetVariable(ds_in, "Ws")
    Wd_in = pfp_utils.GetVariable(ds_in, "Wd")
    # covert to U and V components
    U_in, V_in = pfp_utils.convert_WSWDtoUV(Ws_in, Wd_in)
    # interpolate the components to the output time stamp
    U_out["Data"] = interpolate_1d(x1, U_in["Data"], x2)
    V_out["Data"] = interpolate_1d(x1, V_in["Data"], x2)
    # add the QC flag and update variable attributes
    U_out["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(U_out["Data"])==True, f1, f0)
    V_out["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(V_out["Data"])==True, f1, f0)
    U_out["Attr"]["long_name"] = "U component of wind velocity, postive east"
    U_out["Attr"]["units"] = "m/s"
    V_out["Attr"]["long_name"] = "V component of wind velocity, postive north"
    V_out["Attr"]["units"] = "m/s"
    # write the U and V components to the output data structure
    pfp_utils.CreateVariable(ds_out, U_out)
    pfp_utils.CreateVariable(ds_out, V_out)
    # convert the interpolated components to wind speed and direction
    Ws_out, Wd_out = pfp_utils.convert_UVtoWSWD(U_out, V_out)
    # add the Qc flag and update the variable attributes
    Ws_out["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Ws_out["Data"])==True, f1, f0)
    Wd_out["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Wd_out["Data"])==True, f1, f0)
    Ws_out["Attr"] = copy.deepcopy(Ws_in["Attr"])
    Wd_out["Attr"] = copy.deepcopy(Wd_in["Attr"])
    # write the wind speed and direction into the output data structure
    pfp_utils.CreateVariable(ds_out, Ws_out)
    pfp_utils.CreateVariable(ds_out, Wd_out)
    return

def mask_long_gaps(ds_int, ds_avg, max_interp_length):
    """
    Purpose:
    Assumptions:
    Usage:
    Side effects:
    Author: PRI
    Date: December 2017
    """
    logger.info("Masking long gaps")
    # make a copy of the interpolated data structure
    ds_mlg = copy.deepcopy(ds_int)
    nrecs = ds_int.globalattributes["nc_nrecs"]
    # get a list of the variable labels in this group
    labels = [label for label in ds_avg.series.keys() if label not in ["DateTime","time"]]
    # loop oover the variables in this group
    for label in labels:
        # check to see if this variable exists in the interpolated data structure
        if label not in ds_int.series:
            msg = "mask_long_gaps: variable "+label+" not found, skipping ..."
            logger.warning(msg)
            continue
        # get the average and interpolated variables
        var_avg = pfp_utils.GetVariable(ds_avg, label)
        var_int = pfp_utils.GetVariable(ds_int, label)
        # make a copy of the interpolated variable as a base for the output variable
        var_mlg = copy.deepcopy(var_int)
        # get a boolean array that is True where the average data is masked
        # Note to self: why not use numpy.ma.getmaskarray() here?
        mask = numpy.ma.getmaskarray(var_avg["Data"])
        cond_bool = numpy.ma.where(mask == True, True, False)
        # 2D array of start and end indices for gaps longer that max_interp_length
        cidx = contiguous_regions(cond_bool, max_interp_length = max_interp_length)
        # loop over gaps longer than max_interp_length
        for start, stop in cidx:
            # mask the interpolated data in gaps longer than max_interp_length
            var_mlg["Data"].mask[start:stop] = True
            var_mlg["Flag"][start:stop] = 1
        # write the masked variable back to the output data structure
        pfp_utils.CreateVariable(ds_mlg, var_mlg)
    # return a copy of the interpolated data structure with gaps longer than max_interp_length masked
    return ds_mlg

def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def read_isd_file(isd_file_path):
    """
    Purpose:
     Reads an ISD CSV file (gz or uncompressed) and returns the data in a data structure.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    isd_file_name = os.path.split(isd_file_path)[1]
    msg = "Reading ISD file "+isd_file_name
    logger.info(msg)
    isd_site_id = isd_file_name.split("-")
    isd_site_id = isd_site_id[0]+"-"+isd_site_id[1]
    # read the file
    if os.path.splitext(isd_file_path)[1] == ".gz":
        with gzip.open(isd_file_path, 'rb') as fp:
            content = fp.readlines()
    else:
        with open(isd_file_path) as fp:
            content = fp.readlines()
    # get a data structure
    ds = pfp_io.DataStructure()
    # get the site latitude, longitude and altitude
    ds.globalattributes["altitude"] = float(content[0][46:51])
    ds.globalattributes["latitude"] = float(content[0][28:34])/float(1000)
    ds.globalattributes["longitude"] = float(content[0][34:41])/float(1000)
    ds.globalattributes["isd_site_id"] = isd_site_id
    # initialise the data structure
    isd = {}
    isd["DateTime"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Datetime","units":"none"}}
    isd["Wd"] = {"Data":[],"Attr":{"long_name":"Wind direction","units":"degrees","missing_value":999}}
    isd["Ws"] = {"Data":[],"Attr":{"long_name":"Wind speed","units":"m/s","missing_value":999.9}}
    isd["Ta"] = {"Data":[],"Attr":{"long_name":"Air temperature","units":"C","missing_value":999.9}}
    isd["Td"] = {"Data":[],"Attr":{"long_name":"Dew point temperature","units":"C","missing_value":999.9}}
    isd["ps"] = {"Data":[],"Attr":{"long_name":"Surface pressure","units":"kPa","missing_value":9999.9}}
    isd["Precip"] = {"Data":[],"Attr":{"long_name":"Precipitation","units":"mm","missing_value":999.9}}
    # define the codes for good data in the ISD file
    OK_obs_code = ["AUTO ","CRN05","CRN15","FM-12","FM-15","FM-16","SY-MT"]
    # iterate over the lines in the file and decode the data
    for i in range(len(content)-1):
    #for i in range(10):
        # filter out anything other than hourly data
        if content[i][41:46] not in OK_obs_code: continue
        YY = int(content[i][15:19])
        MM = int(content[i][19:21])
        DD = int(content[i][21:23])
        HH = int(content[i][23:25])
        mm = int(content[i][25:27])
        dt = datetime.datetime(YY,MM,DD,HH,mm,0)
        #isd["DateTime"]["Data"].append(pytz.utc.localize(dt))
        isd["DateTime"]["Data"].append(dt)
        # wind direction, degT
        try:
            isd["Wd"]["Data"].append(float(content[i][60:63]))
        except:
            isd["Wd"]["Data"].append(float(999))
        # wind speed, m/s
        try:
            isd["Ws"]["Data"].append(float(content[i][65:69])/float(10))
        except:
            isd["Ws"]["Data"].append(float(999.9))
        # air temperature, C
        try:
            isd["Ta"]["Data"].append(float(content[i][87:92])/float(10))
        except:
            isd["Ta"]["Data"].append(float(999.9))
        # dew point temperature, C
        try:
            isd["Td"]["Data"].append(float(content[i][93:98])/float(10))
        except:
            isd["Td"]["Data"].append(float(999.9))
        # sea level pressure, hPa
        try:
            isd["ps"]["Data"].append(float(content[i][99:104])/float(10))
        except:
            isd["ps"]["Data"].append(float(9999.9))
        # precipitation, mm
        if content[i][108:111] == "AA1":
            try:
                isd["Precip"]["Data"].append(float(content[i][113:117])/float(10))
            except:
                isd["Precip"]["Data"].append(float(999.9))
        else:
            isd["Precip"]["Data"].append(float(999.9))
    # add the time zone to the DateTime ataributes
    isd["DateTime"]["Attr"]["time_zone"] = "UTC"
    # get the number of records and add this to the global attributes
    nrecs = len(isd["DateTime"]["Data"])
    ds.globalattributes["nc_nrecs"] = str(nrecs)
    # define the QC flags
    f0 = numpy.zeros(len(isd["DateTime"]["Data"]))
    f1 = numpy.ones(len(isd["DateTime"]["Data"]))
    # deal with the datetime first
    variable = {"Label":"DateTime", "Data":numpy.array(isd["DateTime"]["Data"]),
                "Flag":f0, "Attr":isd["DateTime"]["Attr"]}
    pfp_utils.CreateVariable(ds, variable)
    # get the nominal time step
    dt_delta = pfp_utils.get_timestep(ds)
    ts = scipy.stats.mode(dt_delta)[0]/60
    ds.globalattributes["time_step"] = ts[0]
    # add the variables to the data structure
    logger.info("Writing data to the data structure")
    labels = [label for label in isd.keys() if label != "DateTime"]
    for label in labels:
        data = numpy.ma.masked_equal(isd[label]["Data"], isd[label]["Attr"]["missing_value"])
        flag = numpy.where(numpy.ma.getmaskarray(data) == True, f1, f0)
        attr = isd[label]["Attr"]
        variable = {"Label":label, "Data":data, "Flag":flag, "Attr":attr}
        pfp_utils.CreateVariable(ds, variable)
    # hPa to kPa
    ps = pfp_utils.GetVariable(ds, "ps")
    ps["Data"] = ps["Data"]/float(10)
    # convert sea level pressure to station pressure
    site_altitude = float(ds.globalattributes["altitude"])
    Ta = pfp_utils.GetVariable(ds, "Ta")
    cfac = numpy.ma.exp((-1*site_altitude)/((Ta["Data"]+273.15)*29.263))
    ps["Data"] = ps["Data"]*cfac
    ps["Attr"]["long_name"] = ps["Attr"]["long_name"]+", adjusted from sea level to station"
    pfp_utils.CreateVariable(ds, ps)
    # do precipitation and apply crude limits
    Precip = pfp_utils.GetVariable(ds, "Precip")
    condition = (Precip["Data"]<0)|(Precip["Data"]>100)
    Precip["Data"] = numpy.ma.masked_where(condition, Precip["Data"])
    Precip["Flag"] = numpy.where(numpy.ma.getmaskarray(Precip["Data"])==True, f1, f0)
    Precip["Attr"]["RangeCheck_upper"] = 100
    Precip["Attr"]["RangeCheck_lower"] = 0
    pfp_utils.CreateVariable(ds, Precip)
    # get the humidities from Td
    Ta = pfp_utils.GetVariable(ds, "Ta")
    Td = pfp_utils.GetVariable(ds, "Td")
    ps = pfp_utils.GetVariable(ds, "ps")
    RH = mf.RHfromdewpoint(Td["Data"], Ta["Data"])
    flag = numpy.where(numpy.ma.getmaskarray(RH)==True, f1, f0)
    attr = {"long_name":"Relative humidity", "units":"%"}
    variable = {"Label":"RH", "Data":RH, "Flag":flag, "Attr":attr}
    pfp_utils.CreateVariable(ds, variable)
    Ah = mf.absolutehumidityfromRH(Ta["Data"], RH)
    flag = numpy.where(numpy.ma.getmaskarray(Ah)==True, f1, f0)
    attr = {"long_name":"Absolute humidity", "units":"g/m3"}
    variable = {"Label":"Ah", "Data":Ah, "Flag":flag, "Attr":attr}
    pfp_utils.CreateVariable(ds, variable)
    q = mf.specifichumidityfromRH(RH, Ta["Data"], ps["Data"])
    flag = numpy.where(numpy.ma.getmaskarray(q)==True, f1, f0)
    attr = {"long_name":"Specific humidity", "units":"kg/kg"}
    variable = {"Label":"q", "Data":q, "Flag":flag, "Attr":attr}
    pfp_utils.CreateVariable(ds, variable)
    # get U and V components from wind speed and direction
    Ws = pfp_utils.GetVariable(ds, "Ws")
    Wd = pfp_utils.GetVariable(ds, "Wd")
    U, V = pfp_utils.convert_WSWDtoUV(Ws, Wd)
    pfp_utils.CreateVariable(ds, U)
    pfp_utils.CreateVariable(ds, V)
    # add the time variable
    pfp_utils.get_nctime_from_datetime(ds)
    # return the data
    return ds

def read_site_master(xl_file_path, sheet_name):
    """
    """
    xl_book = xlrd.open_workbook(xl_file_path)
    xl_sheet = xl_book.sheet_by_name(sheet_name)
    last_row = int(xl_sheet.nrows)
    # find the header and first data rows
    for i in range(last_row):
        if xl_sheet.cell(i,0).value == "Site":
            header_row = i
            first_data_row = header_row + 1
            break
    # read the header row
    header_row_values = xl_sheet.row_values(header_row)
    # read the site data from the master Excel spreadsheet
    site_info = OrderedDict()
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = OrderedDict()
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value

    return site_info

def toTimestamp(d):
    return calendar.timegm(d.timetuple())

def xl_write_ISD_timesteps(xl_file_path, data):
    """
    Purpose:
     Writes a dictionary to a worksheet in an Excel workbook.
     This routine has 2 arguments,an Excel worksheet instance and
     a dictionary of data to be written out.  The dictionary
     format needs to be:
      data[site][year]["mean"]
      data[site][year]["stdev"]
      data[site][year]["mode"]
    Usage:
     pfp_io.xl_write_ISD_timesteps(xl_file_path, data)
      where xl_file_path is an Excel workbook file name
            data         is a dictionary as defined above
    Side effects:
     Writes to an Excel worksheet instance.
    Called by:
    Calls:
    Author: PRI
    Date: August 2017
    """
    # get a workbook
    xl_book = xlwt.Workbook()
    # get a list of the sheets to add
    site_list = data.keys()
    year_list = data[site_list[0]].keys()
    stat_list = data[site_list[0]][year_list[0]].keys()
    # loop over the statistics
    for stat in stat_list:
        # add a worksheet for the statistics
        xl_sheet = xl_book.add_sheet(stat)
        # write the header line
        for col, year in enumerate(year_list):
            xl_sheet.write(0,col+1,year)
        # write the data, one row per site, one column per year
        for row, site in enumerate(site_list):
            xl_sheet.write(row+1,0,site)
            for col, year in enumerate(year_list):
                if stat in data[site][year].keys():
                    xl_sheet.write(row+1,col+1,data[site][year][stat])
                else:
                    xl_sheet.write(row+1,col+1,"")
    # save the workbook
    xl_book.save(xl_file_path)

    return

# read the control file file
cf = pfp_io.load_controlfile(path='../controlfiles')
xl_file_path = cf["Files"]["xl_file_path"]
xl_sheet_name = cf["Files"]["xl_sheet_name"]
isd_base_path = cf["Files"]["isd_base_path"]
out_base_path = cf["Files"]["out_base_path"]
max_interp_length = cf["Files"]["max_interp_length"]
# read the site master spreadsheet
site_info = read_site_master(xl_file_path, xl_sheet_name)
# get a list of sites
site_list = site_info.keys()
# creat a dictionary to hold the ISD site time steps
isd_time_steps = OrderedDict()

for site in site_list:
    # construct the output file path
    fluxnet_id = site_info[site]["FluxNet ID"]
    if len(fluxnet_id) == 0:
        nc_out_path = os.path.join(out_base_path,site,"Data","ISD",site+"_ISD.nc")
    else:
        nc_out_path = os.path.join(out_base_path,fluxnet_id,"Data","ISD",fluxnet_id+"_ISD.nc")
    # construct the config dictionary for the concatenate routine
    cf_concat = ConfigObj(indent_type="    ")
    cf_concat["Options"] = {"NumberOfDimensions":1,
                            "MaxGapInterpolate":0,
                            "FixTimeStepMethod":"round",
                            "Truncate":"No",
                            "TruncateThreshold":50,
                            "SeriesToCheck":[]}
    cf_concat["Files"] = {"Out":{"ncFileName":nc_out_path}, "In":{}}
    # get the list of ISD stations to be used for this site
    #isd_site_list = cf["Sites"][site]["isd_sites"]
    isd_site_list = []
    for item in ["ISD_ID_1","ISD_ID_2","ISD_ID_3","ISD_ID_4"]:
        if len(site_info[site][item]) != 0:
            isd_site_list.append(site_info[site][item])

    # now get a dictionary that ties the ISD station ID to a number
    # that will be appended to the variable name
    site_index = {}
    for n, isd_site in enumerate(isd_site_list):
        site_index[isd_site] = n
    if not isinstance(isd_site_list, list):
        isd_site_list = [isd_site_list]
    time_zone = site_info[site]["Time zone"]
    time_step = int(round(float(site_info[site]["Time step"])))
    start_year = int(site_info[site]["Start year"])
    end_year = int(site_info[site]["End year"])
    # get the list of years to process
    year_list = range(start_year,end_year+1)
    for n, year in enumerate(year_list):
        # we will collect the data for each site for this year into a single dictionary
        ds_out = {}
        isd_year_path = os.path.join(isd_base_path,str(year))
        for isd_site in isd_site_list:
            if isd_site not in isd_time_steps.keys():
                isd_time_steps[isd_site] = OrderedDict()
            if year not in isd_time_steps[isd_site].keys():
                isd_time_steps[isd_site][year] = OrderedDict()
            isd_file_path = os.path.join(isd_year_path,str(isd_site)+"-"+str(year)+".gz")
            if not os.path.isfile(isd_file_path):
                continue
            ds_in = read_isd_file(isd_file_path)
            pfp_utils.CheckQCFlags(ds_in)
            # get an array of time steps in seconds
            dt = pfp_utils.get_timestep(ds_in)
            # and get dt in minutes
            dt = dt/float(60)
            isd_time_steps[isd_site][year]["mean"] = numpy.mean(dt)
            isd_time_steps[isd_site][year]["stdev"] = numpy.std(dt)
            isd_time_steps[isd_site][year]["mode"] = scipy.stats.mode(dt)[0][0]
            # average records with the same time stamp
            ds_avg = average_duplicate_times(ds_in, time_step)
            pfp_utils.CheckQCFlags(ds_avg)
            # interpolate from the ISD site time step to the tower time step
            ds_int = interpolate_ds(ds_avg, time_step)
            pfp_utils.CheckQCFlags(ds_int)
            # mask long gaps
            ds_mlg = mask_long_gaps(ds_int, ds_avg, max_interp_length)
            pfp_utils.CheckQCFlags(ds_mlg)
            # adjust time from UTC to local using the time zone
            convert_time_zone(ds_mlg, "UTC", time_zone)
            # put the data for this site into the all sites data structure
            ds_out[site_index[isd_site]] = copy.deepcopy(ds_mlg)
            # add some useful global attributes
            ds_out[site_index[isd_site]].globalattributes["isd_site_id"] = isd_site
            ds_out[site_index[isd_site]].globalattributes["time_zone"] = time_zone
            # write out a netCDF file for each ISD site and each year
            #nc_file_name = isd_site+"_"+str(year)+".nc"
            #nc_dir_path = os.path.join(out_base_path,site,"Data","ISD")
            #if not os.path.exists(nc_dir_path):
                #os.makedirs(nc_dir_path)
            #nc_file_path = os.path.join(nc_dir_path,nc_file_name)
            #nc_file = pfp_io.nc_open_write(nc_file_path)
            #pfp_io.nc_write_series(nc_file, ds_out[site_index[isd_site]], ndims=1)
        # now we merge the data structures for each ISD station into a single data structure
        # first, instance a data structure
        ds_all = pfp_io.DataStructure()
        ds_all.globalattributes["latitude"] = site_info[site]["Latitude"]
        ds_all.globalattributes["longitude"] = site_info[site]["Longitude"]
        ds_all.globalattributes["altitude"] = site_info[site]["Altitude"]
        # now loop over the data structures for each ISD station and get the earliest
        # start time and the latest end time
        start_datetime = []
        end_datetime = []
        for i in list(ds_out.keys()):
            start_datetime.append(ds_out[i].series["DateTime"]["Data"][0])
            end_datetime.append(ds_out[i].series["DateTime"]["Data"][-1])
        print site, year
        start = min(start_datetime)
        end = max(end_datetime)
        # now make a datetime series at the required time step from the earliest start
        # datetime to the latest end datetime
        ldt_all = [result for result in perdelta(start, end, datetime.timedelta(minutes=time_step))]
        nrecs = len(ldt_all)
        ds_all.globalattributes["nc_nrecs"] = nrecs
        # and add the datetime to the all-stations data structure
        ds_all.series["DateTime"] = {}
        ds_all.series["DateTime"]["Data"] = ldt_all
        ds_all.series["DateTime"]["Flag"] = numpy.zeros(len(ldt_all))
        ds_all.series["DateTime"]["Attr"] = {"long_name":"DateTime","units":"none"}
        # now copy the contents of the ISD station data structures to the all-stations
        # data structure
        for i in list(ds_out.keys()):
            isd_site_id = ds_out[i].globalattributes["isd_site_id"]
            # copy the global attributes
            gattr_list = list(ds_out[i].globalattributes.keys())
            # remove the site specific global attributes
            for item in ["latitude", "longitude", "altitude", "isd_site_id"]:
                gattr_list.remove(item)
            # copy everything else
            for gattr in gattr_list:
                if gattr not in ds_all.globalattributes:
                    ds_all.globalattributes[gattr] = ds_out[i].globalattributes[gattr]
            # and update the site specific global attributes
            ds_all.globalattributes["time_zone_"+isd_site_id] = ds_out[i].globalattributes["time_zone"]
            ds_all.globalattributes["latitude_"+isd_site_id] = ds_out[i].globalattributes["latitude"]
            ds_all.globalattributes["longitude_"+isd_site_id] = ds_out[i].globalattributes["longitude"]
            ds_all.globalattributes["altitude_"+isd_site_id] = ds_out[i].globalattributes["altitude"]
            # now copy the variables
            # first, we get the indices of matching datetimes
            ldt_one = ds_out[i].series["DateTime"]["Data"]
            idx = numpy.searchsorted(ldt_all, numpy.intersect1d(ldt_all, ldt_one))
            # then we get a list of the variables to copy
            labels = [label for label in ds_out[i].series.keys() if label not in ["DateTime"]]
            # and then we loop over the variables to be copied
            for label in labels:
                # read the data out of the ISD site data structure
                var_out = pfp_utils.GetVariable(ds_out[i], label)
                # create an empty output variable with a number, unique to each ISD station,
                # appended to the label
                var_all = pfp_utils.create_empty_variable(label+"_"+str(i), nrecs)
                # copy the variable attributes
                var_all["Attr"] = copy.deepcopy(var_out["Attr"])
                # add the ISD site ID
                var_all["Attr"]["isd_site_id"] = isd_site_id
                # copy the data and flag onto the matching times
                var_all["Data"][idx] = var_out["Data"]
                var_all["Flag"][idx] = var_out["Flag"]
                # put the data, flag and attributes into the all-in-one data structure
                pfp_utils.CreateVariable(ds_all, var_all)
        # write the netCDF file with the combined data for this year
        if len(fluxnet_id) == 0:
            nc_dir_path = os.path.join(out_base_path,site,"Data","ISD")
            nc_file_name = site+"_ISD_"+str(year)+".nc"
        else:
            nc_dir_path = os.path.join(out_base_path,fluxnet_id,"Data","ISD")
            nc_file_name = fluxnet_id+"_ISD_"+str(year)+".nc"
        if not os.path.exists(nc_dir_path):
            os.makedirs(nc_dir_path)
        nc_file_path = os.path.join(nc_dir_path,nc_file_name)
        nc_file = pfp_io.nc_open_write(nc_file_path)
        pfp_io.nc_write_series(nc_file, ds_all, ndims=1)
        cf_concat["Files"]["In"][str(n)] = nc_file_path
    # concatenate the yearly files for this site
    #cf_concat.filename = "../controlfiles/ISD/concat.txt"
    #cf_concat.write()
    pfp_io.nc_concatenate(cf_concat)

# write the time steps out to an Excel file
xl_file_path = os.path.join(isd_base_path, "ISD_site_timesteps.xls")
xl_write_ISD_timesteps(xl_file_path, isd_time_steps)

logger.info("All done")
