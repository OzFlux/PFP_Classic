# standard modules
import copy
import os
import sys
# 3rd party modules
from configobj import ConfigObj
import timezonefinder
# check the scripts folder exists
scripts_path = os.path.join("..", "scripts", "")
if not os.path.exists(scripts_path):
    print "cleanup_netcdf_files: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append(scripts_path)
# PFP modules
import constants as c
import pfp_io
import pfp_utils

def change_variable_names(cfg, ds):
    """
    Purpose:
     Change variable names to the new (October 2018) scheme.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # get a list of potential mappings
    rename_list = [v for v in list(cfg["Variables"].keys()) if "rename" in cfg["Variables"][v]]
    # loop over the variables in the data structure
    series_list = list(ds.series.keys())
    for label in series_list:
        if label in rename_list:
            new_name = cfg["Variables"][label]["rename"]
            ds.series[new_name] = ds.series.pop(label)
    return

def copy_ws_wd(ds):
    """
    Purpose:
     Make sure the Ws and Wd variables are in the L3 netCDF files.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # get a list of the series
    series_list = sorted(list(ds.series.keys()))
    if "Wd" not in series_list:
        if "Wd_SONIC_Av" in series_list:
            ds.series["Wd"] = copy.deepcopy(ds.series["Wd_SONIC_Av"])
            ds.series["Wd"]["Attr"]["long_name"] = "Wind direction (copied from Wd_SONIC_Av)"
    if "Ws" not in series_list:
        if "Ws_SONIC_Av" in series_list:
            ds.series["Ws"] = copy.deepcopy(ds.series["Ws_SONIC_Av"])
            ds.series["Ws"]["Attr"]["long_name"] = "Wind speed (copied from Ws_SONIC_Av)"
    return

def remove_variables(cfg, ds):
    """
    Purpose:
     Remove deprecated variables from a netCDF file.
    Usage:
    Author: PRI
    Date: October 2018
    """
    remove_list = [v for v in list(cfg["Variables"].keys()) if "remove" in cfg["Variables"][v]]
    series_list = sorted(list(ds.series.keys()))
    for label in series_list:
        if label in remove_list:
            ds.series.pop(label)
    return

def change_global_attributes(cfg, ds):
    """
    Purpose:
     Clean up the global attributes.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # check site_name is in ds.globalattributes
    gattr_list = list(ds.globalattributes.keys())
    if "site_name" not in gattr_list:
        print "Globel attributes: site_name not found"
    # check latitude and longitude are in ds.globalattributes
    if "latitude" not in gattr_list:
        print "Global attributes: latitude not found"
    else:
        lat_string = str(ds.globalattributes["latitude"])
        if len(lat_string) == 0:
            print "Global attributes: latitude empty"
        else:
            lat = pfp_utils.convert_anglestring(lat_string)
        ds.globalattributes["latitude"] = str(lat)
    if "longitude" not in gattr_list:
        print "Global attributes: longitude not found"
    else:
        lon_string = str(ds.globalattributes["longitude"])
        if len(lon_string) == 0:
            print "Global attributes: longitude empty"
        else:
            lon = pfp_utils.convert_anglestring(lon_string)
        ds.globalattributes["longitude"] = str(lon)
    # check to see if there there is a time_zone global attribute
    gattr_list = list(ds.globalattributes.keys())
    if not "time_zone" in gattr_list:
        # get the site name
        site_name = ds.globalattributes["site_name"]
        sn = site_name.replace(" ","").replace(",","").lower()
        # first, see if the site is in constants.tz_dict
        if sn in list(c.tz_dict.keys()):
            ds.globalattributes["time_zone"] = c.tz_dict[sn]
        else:
            if "latitude" in gattr_list and "longitude" in gattr_list:
                lat = float(ds.globalattributes["latitude"])
                lon = float(ds.globalattributes["longitude"])
                if lat != -9999 and lon != -9999:
                    tf = timezonefinder.TimezoneFinder()
                    tz = tf.timezone_at(lng=lon, lat=lat)
                    ds.globalattributes["time_zone"] = tz
                else:
                    print "Global attributes: unable to define time zone"
                    ds.globalattributes["time_zone"] = ""
    # add or change global attributes as required
    gattr_list = sorted(list(cfg["Global"].keys()))
    for gattr in gattr_list:
        ds.globalattributes[gattr] = cfg["Global"][gattr]
    # remove deprecated global attributes
    flag_list = [g for g in ds.globalattributes.keys() if "Flag" in g]
    others_list = ["end_datetime", "start_datetime", "Functions", "doi"]
    remove_list = others_list + flag_list
    for gattr in list(ds.globalattributes.keys()):
        if gattr in remove_list:
            ds.globalattributes.pop(gattr)
    # rename global attributes
    rename_dict = {"EPDversion":"PythonVersion", "elevation":"altitude"}
    for item in rename_dict:
        if item in list(ds.globalattributes.keys()):
            new_key = rename_dict[item]
            ds.globalattributes[new_key] = ds.globalattributes.pop(item)
    return

cfg_name = os.path.join("..", "controlfiles", "standard", "map_old_to_new.txt")
if os.path.exists(cfg_name):
    cfg = ConfigObj(cfg_name)
else:
    print " 'map_old_to_new' control file not found"

rp = os.path.join(os.sep, "mnt", "OzFlux", "Sites")
sites = sorted([d for d in os.listdir(rp) if os.path.isdir(os.path.join(rp,d))])

for site in sites:
    sp = os.path.join(rp, site, "Data", "Portal")
    if not os.path.isdir(sp):
        print sp+" , skipping site ..."
        continue
    files = sorted([f for f in os.listdir(sp) if ("L3" in f and ".nc" in f)])
    if len(files) == 0:
        print "No files found in "+sp+" , skipping ..."
        continue
    for file in files:
        fp = os.path.join(sp, file)
        print "Converting "+file
        ds = pfp_io.nc_read_series(fp)
        change_variable_names(cfg, ds)
        copy_ws_wd(ds)
        remove_variables(cfg, ds)
        change_global_attributes(cfg, ds)
        nf = pfp_io.nc_open_write(fp)
        pfp_io.nc_write_series(nf, ds)