import ast
import datetime
import logging
import ntpath
import os
import sys
sys.path.append('scripts')
import time
import pfp_clim
import pfp_cpd
import pfp_io
import pfp_log
import pfp_ls
import pfp_plot
import pfp_utils

t = time.localtime()
rundatetime = datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]).strftime("%Y%m%d%H%M")
log_filename = 'batchprocess_'+rundatetime+'.log'

logger = pfp_log.init_logger(logger_name="pfp_log", file_handler=log_filename)

# get the batch processing control file
if len(sys.argv)==1:
    cf_batch = pfp_io.load_controlfile(path='controlfiles')
    if len(cf_batch)==0: sys.exit()
else:
    cfname = sys.argv[1]
    if os.path.exists(cfname):
        cf_batch = pfp_io.get_controlfilecontents(cfname)
    else:
        logger.error("Control file "+cfname+" does not exist")
        sys.exit()

level_list = ['L1','L2','L3','concatenate','climatology','cpd','L4','L5','L6']
if "Options" in cf_batch:
    if "levels" in cf_batch["Options"]: level_list = ast.literal_eval(cf_batch["Options"]["levels"])
for level in level_list:
    if level.lower()=="l1":
        # L1 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L1 processing with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            ds1 = pfp_ls.l1qc(cf)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds1)
            logger.info('Finished L1 processing with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="l2":
        # L2 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L2 processing with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            infilename = pfp_io.get_infilenamefromcf(cf)
            ds1 = pfp_io.nc_read_series(infilename)
            ds2 = pfp_ls.l2qc(cf,ds1)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds2)
            logger.info('Finished L2 processing with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="l3":
        # L3 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L3 processing with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            infilename = pfp_io.get_infilenamefromcf(cf)
            ds2 = pfp_io.nc_read_series(infilename)
            ds3 = pfp_ls.l3qc(cf,ds2)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            outputlist = pfp_io.get_outputlistfromcf(cf,'nc')
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds3,outputlist=outputlist)
            logger.info('Finished L3 processing with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="fluxnet":
        # convert netCDF files to FluxNet CSV files
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            cf_file_name = os.path.split(cfname)
            logger.info('Starting FluxNet output with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            pfp_io.fn_write_csv(cf)
            logger.info('Finished FluxNet output with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="reddyproc":
        # convert netCDF files to REddyProc CSV files
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            cf_file_name = os.path.split(cfname)
            logger.info('Starting REddyProc output with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            pfp_io.reddyproc_write_csv(cf)
            logger.info('Finished REddyProc output with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="concatenate":
        # concatenate netCDF files
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting concatenation with '+cf_file_name[1])
            cf_cc = pfp_io.get_controlfilecontents(cfname)
            pfp_io.nc_concatenate(cf_cc)
            logger.info('Finished concatenation with '+cf_file_name[1])
            # now plot the fingerprints for the concatenated files
            opt = pfp_utils.get_keyvaluefromcf(cf_cc,["Options"],"DoFingerprints", default="yes")
            if opt.lower()=="no": continue
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp): cf_fp["Files"] = {}
            file_name = cf_cc["Files"]["Out"]["ncFileName"]
            file_path = ntpath.split(file_name)[0]+"/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")]+"Plots/"
            if "Options" not in cf_fp: cf_fp["Options"]={}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "no"
            logger.info('Doing fingerprint plots using '+cf_fp["Files"]["in_filename"])
            pfp_plot.plot_fingerprint(cf_fp)
            logger.info('Finished fingerprint plots')
            logger.info('')
    elif level.lower()=="climatology":
        # climatology
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting climatology with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            pfp_clim.climatology(cf)
            logger.info('Finished climatology with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="cpd":
        # ustar threshold from change point detection
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting CPD with '+cf_file_name[1])
            cf = pfp_io.get_controlfilecontents(cfname)
            if "Options" not in cf: cf["Options"]={}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = False
            pfp_cpd.cpd_main(cf)
            logger.info('Finished CPD with '+cf_file_name[1])
            logger.info('')
    elif level.lower()=="l4":
        # L4 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L4 processing with '+cf_file_name[1])
            cf_l4 = pfp_io.get_controlfilecontents(cfname)
            if "Options" not in cf_l4: cf_l4["Options"]={}
            cf_l4["Options"]["call_mode"] = "batch"
            cf_l4["Options"]["show_plots"] = False
            infilename = pfp_io.get_infilenamefromcf(cf_l4)
            ds3 = pfp_io.nc_read_series(infilename)
            ds4 = pfp_ls.l4qc(cf_l4,ds3)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l4)
            outputlist = pfp_io.get_outputlistfromcf(cf_l4,'nc')
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds4,outputlist=outputlist)
            logger.info('Finished L4 processing with '+cf_file_name[1])
            # now plot the fingerprints for the L4 files
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp): cf_fp["Files"] = {}
            file_name = pfp_io.get_outfilenamefromcf(cf_l4)
            file_path = ntpath.split(file_name)[0]+"/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")]+"Plots/"
            if "Options" not in cf_fp: cf_fp["Options"]={}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "no"
            logger.info('Doing fingerprint plots using '+cf_fp["Files"]["in_filename"])
            pfp_plot.plot_fingerprint(cf_fp)
            logger.info('Finished fingerprint plots')
            logger.info('')
    elif level.lower()=="l5":
        # L5 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L5 processing with '+cf_file_name[1])
            cf_l5 = pfp_io.get_controlfilecontents(cfname)
            if "Options" not in cf_l5: cf_l5["Options"]={}
            cf_l5["Options"]["call_mode"] = "batch"
            cf_l5["Options"]["show_plots"] = False
            infilename = pfp_io.get_infilenamefromcf(cf_l5)
            ds4 = pfp_io.nc_read_series(infilename)
            ds5 = pfp_ls.l5qc(cf_l5,ds4)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l5)
            outputlist = pfp_io.get_outputlistfromcf(cf_l5,'nc')
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds5,outputlist=outputlist)
            logger.info('Finished L5 processing with '+cf_file_name[1])
            # now plot the fingerprints for the L5 files
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp): cf_fp["Files"] = {}
            file_name = pfp_io.get_outfilenamefromcf(cf_l5)
            file_path = ntpath.split(file_name)[0]+"/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")]+"Plots/"
            if "Options" not in cf_fp: cf_fp["Options"]={}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "no"
            logger.info('Doing fingerprint plots using '+cf_fp["Files"]["in_filename"])
            pfp_plot.plot_fingerprint(cf_fp)
            logger.info('Finished fingerprint plots')
            logger.info('')
    elif level.lower()=="l6":
        # L6 processing
        for i in cf_batch["Levels"][level].keys():
            cfname = cf_batch["Levels"][level][i]
            if not os.path.isfile(cfname):
                msg = " Control file "+cfname+" not found"
                logger.error(msg)
                continue
            cf_file_name = os.path.split(cfname)
            logger.info('Starting L6 processing with '+cf_file_name[1])
            try:
                cf = pfp_io.get_controlfilecontents(cfname)
                if "Options" not in cf: cf["Options"]={}
                cf["Options"]["call_mode"] = "batch"
                cf["Options"]["show_plots"] = False
                infilename = pfp_io.get_infilenamefromcf(cf)
                ds5 = pfp_io.nc_read_series(infilename)
                ds6 = pfp_ls.l6qc(cf,ds5)
                outfilename = pfp_io.get_outfilenamefromcf(cf)
                outputlist = pfp_io.get_outputlistfromcf(cf,'nc')
                ncFile = pfp_io.nc_open_write(outfilename)
                pfp_io.nc_write_series(ncFile,ds6,outputlist=outputlist)
                logger.info('Finished L6 processing with '+cf_file_name[1])
                logger.info('')
            except:
                msg = "Error occurred during L6 with "+cf_file_name[1]
                logger.error(msg)
                continue
