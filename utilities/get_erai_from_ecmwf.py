#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import datetime
import dateutil.parser
import sys

erai_info = {}
erai_info["class"] = "ei"
erai_info["dataset"] = "interim"
erai_info["expver"] = "1"
erai_info["grid"] = "0.75/0.75"
erai_info["levtype"] = "sfc"
erai_info["param"] = "39.128/40.128/41.128/42.128/134.128/139.128/146.128/147.128/159.128/165.128/166.128/167.128/168.128/169.128/170.128/175.128/176.128/177.128/183.128/228.128/236.128"
erai_info["step"] = "3/6/9/12"
erai_info["stream"] = "oper"
erai_info["time"] = "00/12"
erai_info["type"] = "fc"
erai_info["format"] = "netcdf"

if len(sys.argv)==1:
    print "Command line syntax is:"
    print " python get_erai.py <country>"
    print "where <country> can be;"
    print " Australia"
    print " USA"
    print " NZ"
    sys.exit

if sys.argv[1].lower()=="australia":
    print "Get ERAI data for Australia"
    erai_info["area"] = "-10/110/-45/155"
    target_directory = "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/AUS/"
    start_date = "2018-01-01"
    end_date = "2018-04-30"
    print "Start=",start_date," End=",end_date
elif sys.argv[1].lower()=="nz":
    print "Get ERAI data for New Zealand"
    erai_info["area"] = "-30/165/-50/180"
    target_directory = "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/NZ/"
    start_date = "2011-01-01"
    end_date = "2018-02-28"
    print "Start=",start_date," End=",end_date
elif sys.argv[1].lower()=="usa":
    erai_info["area"] = "70/229.5/30/300"
    target_directory = "/home/peter/AmeriFlux/ERAI/"
    start_date = "2017-01-01"
    end_date = "2017-10-31"
elif sys.argv[1].lower()=="china":
    erai_info["area"] = "60/80/20/140"
    target_directory = "/home/peter/ChinaFlux/ERAI/"
    start_date = "2014-01-01"
    end_date = "2014-12-31"
elif sys.argv[1].lower()=="borneo":
    print "Get ERAI data for Borneo/Indonesia"
    erai_info["area"] = "8.25/108/3.75/120"
    target_directory = "/run/media/cilli/cillidata/cilli/1_Work/1_OzFlux/Sites/ERAI/BORNEO/"
    start_date = "2011-01-01"
    end_date = "2018-05-31"
    print "Start=",start_date," End=",end_date
else:
    print "Unrecognised country option entered on command line"
    print "Valid country options are:"
    print " australia"
    print " usa"
    sys.exit()

server = ECMWFDataServer()
sd = dateutil.parser.parse(start_date)
ed = dateutil.parser.parse(end_date)
start_year = sd.year
end_year = ed.year
year_list = range(start_year,end_year+1)
for year in year_list:
    print " Processing year: ",str(year)
    sds = str(year)+"-01-01"
    edc = datetime.datetime(year+1,1,1,0,0,0)
    eds = edc.strftime("%Y-%m-%d")
    if ed < edc:
        eds = ed.strftime("%Y-%m-%d")
    erai_info["date"] = sds+"/to/"+eds
    #print sds+"/to/"+eds
    erai_info["target"] = target_directory+"ERAI_"+str(year)+".nc"
    server.retrieve(erai_info)
