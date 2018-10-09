#import ast
#import copy
import datetime
#import logging
import matplotlib
#import matplotlib.style
matplotlib.use('TkAgg')
#matplotlib.style.use('classic')
#matplotlib.use('Qt4Agg')
#import numpy
import ntpath
import time
import Tkinter as tk
import tkMessageBox
import os
import sys

# The Lindsay Trap: check the scripts directory is present
if not os.path.exists("./scripts/"):
    print "OzFluxQC: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('scripts')
import cfg
import pfp_clim
import pfp_cpd
import pfp_gf
import pfp_io
import pfp_log
import pfp_ls
import pfp_mpt
import pfp_plot
import pfp_rp
import pfp_utils
# now check the logfiles and plots directories are present
dir_list = ["./logfiles/", "./plots/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# now check the solo/inf, solo/input, solo/log and solo/output directories are present
dir_list = ["./solo/inf", "./solo/input", "./solo/log", "./solo/output"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# next we make sure the MPT directories are present ...
dir_list = ["./mpt/input", "./mpt/log", "./mpt/output"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# ... and make sure the MDS directories are present
dir_list = ["./mds/input", "./mds/log", "./mds/output"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# start a log file with the current date and time in the name
t = time.localtime()
rundatetime = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5]).strftime("%Y%m%d%H%M")
log_filename = 'pfp_'+rundatetime+'.log'
logger = pfp_log.init_logger(logger_name="pfp_log", file_handler=log_filename)

class qcgui(tk.Tk):
    """
        QC Data Main GUI
        Used to access read, save, and data processing (pfp_ls) prodecures

        Columns: Data levels:
            1:  L1 Raw Data (read excel into NetCDF)
            2:  L2 QA/QC (general QA/QC algorithms, site independent)
            3:  L3 Corrections (Flux data corrections, site dependent based on ancillary measurements available and technical issues)
            4:  L4 Gap Filling (Used for fill met data gaps and ingesting SOLO-ANN Gap Filled fluxes from external processes)

        Rows:  function access
            1:  Ingest excel dataset into NetCDF files
            2:  Process data from previous level and generate NetCDF file(s) at current level
            3-6:  Show Timestamp range of dataset and accept date range for graphical plots
            7:  Export excel dataset from NetCDF file
        """
    def __init__(self, parent):
        tk.Tk.__init__(self, parent)
        self.parent = parent
        self.initialise()

    def option_not_implemented(self):
        self.do_progress(text='Option not implemented yet ...')
        logger.info(' Option not implemented yet ...')

    def initialise(self):
        self.org_frame = tk.Frame(self)
        self.org_frame.grid()
        # things in the first row of the GUI
        L1Label = tk.Label(self.org_frame, text='L1: Raw data')
        L1Label.grid(row=0, column=0, columnspan=2)
        L2Label = tk.Label(self.org_frame, text='L2: QA/QC')
        L2Label.grid(row=0, column=2, columnspan=2)
        L3Label = tk.Label(self.org_frame, text='L3: Process')
        L3Label.grid(row=0, column=4, columnspan=2)
        # things in the second row of the GUI
        doL1Button = tk.Button(self.org_frame, text="Read L1 file", command=self.do_l1qc)
        doL1Button.grid(row=1, column=0, columnspan=2)
        doL2Button = tk.Button(self.org_frame, text="Do L2 QA/QC", command=self.do_l2qc)
        doL2Button.grid(row=1, column=2, columnspan=2)
        doL3Button = tk.Button(self.org_frame, text="Do L3 processing", command=self.do_l3qc)
        doL3Button.grid(row=1, column=4, columnspan=2)
        # things in the third row of the GUI
        filestartLabel = tk.Label(self.org_frame,text='File start date')
        filestartLabel.grid(row=2,column=0,columnspan=3)
        fileendLabel = tk.Label(self.org_frame,text='File end date')
        fileendLabel.grid(row=2,column=3,columnspan=3)
        # things in the fourth row of the GUI
        self.filestartValue = tk.Label(self.org_frame,text='No file loaded ...')
        self.filestartValue.grid(row=3,column=0,columnspan=3)
        self.fileendValue = tk.Label(self.org_frame,text='No file loaded ...')
        self.fileendValue.grid(row=3,column=3,columnspan=3)
        # things in the fifth row of the GUI
        plotstartLabel = tk.Label(self.org_frame, text='Start date (YYYY-MM-DD)')
        plotstartLabel.grid(row=4,column=0,columnspan=3)
        self.plotstartEntry = tk.Entry(self.org_frame)
        self.plotstartEntry.grid(row=4,column=3,columnspan=3)
        # things in row sixth of the GUI
        plotendLabel = tk.Label(self.org_frame, text='End date   (YYYY-MM-DD)')
        plotendLabel.grid(row=5,column=0,columnspan=3)
        self.plotendEntry = tk.Entry(self.org_frame)
        self.plotendEntry.grid(row=5,column=3,columnspan=3)
        # things in the seventh row of the GUI
        closeplotwindowsButton = tk.Button (self.org_frame, text="Close plot windows", command=self.do_closeplotwindows )
        closeplotwindowsButton.grid(row=6,column=0,columnspan=2)
        plotL1L2Button = tk.Button (self.org_frame, text="Plot L1 & L2 Data", command=self.do_plotL1L2 )
        plotL1L2Button.grid(row=6,column=2,columnspan=2)
        plotL3L3Button = tk.Button (self.org_frame, text="Plot L3 Data", command=self.do_plotL3L3 )
        plotL3L3Button.grid(row=6,column=4,columnspan=2)
        # things in the eigth row of the GUI
        quitButton = tk.Button (self.org_frame, text='Quit', command=self.do_quit )
        quitButton.grid(row=7,column=0,columnspan=2)
        savexL2Button = tk.Button (self.org_frame, text='Write L2 Excel file', command=self.do_savexL2 )
        savexL2Button.grid(row=7,column=2,columnspan=2)
        savexL3Button = tk.Button (self.org_frame, text='Write L3 Excel file', command=self.do_savexL3 )
        savexL3Button.grid(row=7,column=4,columnspan=2)
        # other things in the GUI
        self.progress = tk.Label(self.org_frame, text='Waiting for input ...')
        self.progress.grid(row=8,column=0,columnspan=6,sticky="W")
        # now we put together the menu, "File" first
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar,tearoff=0)
        filemenu.add_command(label="Concatenate netCDF",command=self.do_ncconcat)
        filemenu.add_command(label="Split netCDF",command=self.do_ncsplit)
        filemenu.add_command(label="List netCDF contents",command=self.option_not_implemented)
        fileconvertmenu = tk.Menu(menubar,tearoff=0)
        #fileconvertmenu.add_command(label="V2.7 to V2.8",command=self.do_v27tov28)
        fileconvertmenu.add_command(label="nc to EddyPro (biomet)",command=self.do_nc2ep_biomet)
        fileconvertmenu.add_command(label="nc to ECOSTRESS",command=self.do_nc2csv_ecostress)
        fileconvertmenu.add_command(label="nc to FluxNet",command=self.do_nc2fn)
        fileconvertmenu.add_command(label="nc to REddyProc",command=self.do_nc2reddyproc)
        fileconvertmenu.add_command(label="nc to SMAP",command=self.do_nc2smap)
        fileconvertmenu.add_command(label="nc to xls",command=self.do_nc2xls)
        fileconvertmenu.add_command(label="xls to nc",command=self.option_not_implemented)
        filemenu.add_cascade(label="Convert",menu=fileconvertmenu)
        filemenu.add_separator()
        filemenu.add_command(label="Quit",command=self.do_quit)
        menubar.add_cascade(label="File",menu=filemenu)
        # now the "Run" menu
        runmenu = tk.Menu(menubar,tearoff=0)
        runmenu.add_command(label="Read L1 file",command=self.do_l1qc)
        runmenu.add_command(label="Do L2 QA/QC",command=self.do_l2qc)
        runmenu.add_command(label="Do L3 processing",command=self.do_l3qc)
        runmenu.add_command(label="Do L4 gap fill (drivers)",command=self.do_l4qc)
        runmenu.add_command(label="Do L5 gap fill (fluxes)",command=self.do_l5qc)
        runmenu.add_command(label="Do L6 partitioning",command=self.do_l6qc)
        #runmenu.add_command(label="Batch processing",command=self.do_batch)
        menubar.add_cascade(label="Run",menu=runmenu)
        # then the "Plot" menu
        plotmenu = tk.Menu(menubar,tearoff=0)
        plotmenu.add_command(label="Plot L1 & L2",command=self.do_plotL1L2)
        plotmenu.add_command(label="Plot L3",command=self.do_plotL3L3)
        plotmenu.add_command(label="Plot L4",command=self.do_plotL3L4)
        plotmenu.add_command(label="Plot L5",command=self.option_not_implemented)
        plotmenu.add_command(label="Plot L6 summary",command=self.do_plotL6_summary)
        fnmenu = tk.Menu(menubar,tearoff=0)
        fnmenu.add_command(label="Standard",command=lambda:self.do_plotfluxnet(mode="standard"))
        fnmenu.add_command(label="Custom",command=lambda:self.do_plotfluxnet(mode="custom"))
        plotmenu.add_cascade(label="30 minute",menu=fnmenu)
        #plotmenu.add_command(label="FluxNet",command=self.do_plotfluxnet)
        fpmenu = tk.Menu(menubar,tearoff=0)
        fpmenu.add_command(label="Standard",command=lambda:self.do_plotfingerprint(mode="standard"))
        fpmenu.add_command(label="Custom",command=lambda:self.do_plotfingerprint(mode="custom"))
        plotmenu.add_cascade(label="Fingerprint",menu=fpmenu)
        plotmenu.add_command(label="Quick check",command=self.do_plotquickcheck)
        plotmenu.add_command(label="Years check",command=self.option_not_implemented)
        plotmenu.add_separator()
        plotmenu.add_command(label="Close plots",command=self.do_closeplotwindows)
        menubar.add_cascade(label="Plot",menu=plotmenu)
        # and the "Utilities" menu
        utilsmenu = tk.Menu(menubar,tearoff=0)
        # climatology menu
        climatologymenu = tk.Menu(menubar,tearoff=0)
        climatologymenu.add_command(label="Standard",command=lambda:self.do_climatology(mode="standard"))
        climatologymenu.add_command(label="Custom",command=lambda:self.do_climatology(mode="custom"))
        utilsmenu.add_cascade(label="Climatology",menu=climatologymenu)
        # CPD menu
        cpdmenu = tk.Menu(menubar,tearoff=0)
        cpdmenu.add_command(label="Standard",command=lambda:self.do_cpd(mode="standard"))
        cpdmenu.add_command(label="Custom",command=lambda:self.do_cpd(mode="custom"))
        utilsmenu.add_cascade(label="u* threshold (CPD)",menu=cpdmenu)
        mptmenu = tk.Menu(menubar,tearoff=0)
        mptmenu.add_command(label="Standard",command=lambda:self.do_mpt(mode="standard"))
        mptmenu.add_command(label="Custom",command=lambda:self.do_mpt(mode="custom"))
        utilsmenu.add_cascade(label="u* threshold (MPT)",menu=mptmenu)
        menubar.add_cascade(label="Utilities",menu=utilsmenu)
        # and the "Help" menu
        helpmenu = tk.Menu(menubar,tearoff=0)
        helpmenu.add_command(label="Contents",command=self.do_helpcontents)
        helpmenu.add_command(label="About",command=self.option_not_implemented)
        menubar.add_cascade(label="Help",menu=helpmenu)

        self.config(menu=menubar)

    def do_climatology(self,mode="standard"):
        """
        Calls pfp_clim.climatology
        """
        logger.info(' Starting climatology')
        self.do_progress(text='Doing climatology ...')
        if mode=="standard":
            stdname = "controlfiles/standard/climatology.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                self.do_progress(text='Opening input file ...')
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if not os.path.exists(filename):
                    logger.info( " Climatology: no input file chosen")
                    self.do_progress(text='Waiting for input ...')
                    return
                if "Files" not in dir(cf): cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                in_filename = ntpath.split(filename)[1]
                cf["Files"]["in_filename"] = in_filename
                cf["Files"]["out_filename"] = in_filename.replace(".nc","_Climatology.xls")
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Doing the climatology')
        pfp_clim.climatology(cf)
        self.do_progress(text='Finished climatology')
        logger.info(' Finished climatology')
        logger.info("")

    def do_closeplotwindows(self):
        """
            Close plot windows
            """
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        logger.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        #fig_numbers = [n.num for n in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        ##logger.info('  Closing plot windows: '+str(fig_numbers))
        #for n in fig_numbers:
            #matplotlib.pyplot.close(n)
        self.do_progress(text='Waiting for input ...')             # tell the user what we're doing
        logger.info(' Waiting for input ...')

    def do_compare_eddypro(self):
        """
        Calls pfp_clim.compare_ep
        Compares the results OzFluxQC (L3) with those from EddyPro (full output).
        """
        self.do_progress(text='Comparing EddyPro and OzFlux results ...')
        pfp_clim.compare_eddypro()
        self.do_progress(text='Finished comparing EddyPro and OzFlux')
        logger.info(' Finished comparing EddyPro and OzFlux')

    def do_cpd(self,mode="standard"):
        """
        Calls pfp_cpd.cpd_main
        Calculate the u* threshold using the Change Point Detection (CPD) method.
        """
        logger.info(' Starting estimation u* threshold using CPD')
        self.do_progress(text='Estimating u* threshold using CPD ...')
        if mode=="standard":
            stdname = "controlfiles/standard/cpd.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                self.do_progress(text='Opening input file ...')
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if not os.path.exists(filename):
                    logger.info( " CPD: no input file chosen")
                    self.do_progress(text='Waiting for input ...')
                    return
                if "Files" not in dir(cf): cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                in_filename = ntpath.split(filename)[1]
                cf["Files"]["in_filename"] = in_filename
                cf["Files"]["out_filename"] = in_filename.replace(".nc","_CPD.xls")
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf)==0:
                    self.do_progress(text='Waiting for input ...')
                    return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        self.do_progress(text='Doing the u* threshold (CPD)')
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        pfp_cpd.cpd_main(cf)
        self.do_progress(text='Finished estimating u* threshold')
        logger.info(' Finished estimating u* threshold')
        logger.info("")

    def do_helpcontents(self):
        tkMessageBox.showinfo("Obi Wan says ...","Read the source, Luke!")

    def do_l1qc(self):
        """
        Calls pfp_ls.l1qc
        """
        logger.info(" Starting L1 processing ...")
        self.do_progress(text='Load L1 Control File ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            logger.info( " L1: no control file chosen")
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Doing L1 ...')
        ds1 = pfp_ls.l1qc(cf)
        if ds1.returncodes["value"] == 0:
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            ncFile = pfp_io.nc_open_write(outfilename)
            pfp_io.nc_write_series(ncFile,ds1)
            self.do_progress(text='Finished L1')
            logger.info(' Finished L1')
            logger.info("")
        else:
            msg = 'An error occurred, check the console ...'
            self.do_progress(text=msg)

    def do_l2qc(self):
        """
            Call pfp_ls.l2qc function
            Performs L2 QA/QC processing on raw data
            Outputs L2 netCDF file to ncData folder

            ControlFiles:
                L2_year.txt
                or
                L2.txt

            ControlFile contents (see ControlFile/Templates/L2.txt for example):
                [General]:
                    Enter list of functions to be performed
                [Files]:
                    L1 input file name and path
                    L2 output file name and path
                [Variables]:
                    Variable names and parameters for:
                        Range check to set upper and lower rejection limits
                        Diurnal check to reject observations by time of day that
                            are outside specified standard deviation limits
                        Timestamps for excluded dates
                        Timestamps for excluded hours
                [Plots]:
                    Variable lists for plot generation
            """
        logger.info(" Starting L2 processing ...")
        self.do_progress(text='Load L2 Control File ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            logger.info( " L2: no control file chosen")
            self.do_progress(text='Waiting for input ...')
            return
        infilename = pfp_io.get_infilenamefromcf(cf)
        if not pfp_utils.file_exists(infilename):
            self.do_progress(text='An error occurred, check the console ...')
            return
        self.do_progress(text='Doing L2 QC ...')
        self.ds1 = pfp_io.nc_read_series(infilename)
        if len(self.ds1.series.keys())==0:
            self.do_progress(text='An error occurred, check the console ...')
            del self.ds1
            return
        self.update_startenddate(str(self.ds1.series['DateTime']['Data'][0]),
                                 str(self.ds1.series['DateTime']['Data'][-1]))
        self.ds2 = pfp_ls.l2qc(cf, self.ds1)
        # put a copy of the config file in self for later use
        self.cf = cf
        logger.info(' Finished L2 QC process')
        self.do_progress(text='Finished L2 QC process')
        self.do_progress(text='Saving L2 QC ...')                     # put up the progress message
        outfilename = pfp_io.get_outfilenamefromcf(cf)
        if len(outfilename)==0:
            self.do_progress(text='An error occurred, check the console ...')
            return
        ncFile = pfp_io.nc_open_write(outfilename)
        pfp_io.nc_write_series(ncFile, self.ds2)                                  # save the L2 data
        self.do_progress(text='Finished saving L2 QC data')              # tdo_progressell the user we are done
        logger.info(' Finished saving L2 QC data')
        logger.info("")

    def do_l3qc(self):
        """
            Call pfp_ls.l3qc_sitename function
            Performs L3 Corrections and QA/QC processing on L2 data
            Outputs L3 netCDF file to ncData folder
            Outputs L3 netCDF file to OzFlux folder

            Available corrections:
            * corrections requiring ancillary measurements or samples
              marked with an asterisk
                Linear correction
                    fixed slope
                    linearly shifting slope
                Conversion of virtual temperature to actual temperature
                2D Coordinate rotation
                Massman correction for frequency attenuation*
                Webb, Pearman and Leuning correction for flux effects on density
                    measurements
                Conversion of virtual heat flux to actual heat flux
                Correction of soil moisture content to empirical calibration
                    curve*
                Addition of soil heat storage to ground ground heat flux*

            ControlFiles:
                L3_year.txt
                or
                L3a.txt

            ControlFile contents (see ControlFile/Templates/L3.txt for example):
                [General]:
                    Python control parameters
                [Files]:
                    L2 input file name and path
                    L3 output file name and ncData folder path
                    L3 OzFlux output file name and OzFlux folder path
                [Massman] (where available):
                    Constants used in frequency attenuation correction
                        zmd: instrument height (z) less zero-plane displacement
                            height (d), m
                        z0: aerodynamic roughness length, m
                        angle: angle from CSAT mounting point between CSAT and
                            IRGA mid-path, degrees
                        CSATarm: distance from CSAT mounting point to CSAT
                            mid-path, m
                        IRGAarm: distance from CSAT mounting point to IRGA
                            mid-path, m
                [Soil]:
                    Constants used in correcting Fg for storage and in empirical
                    corrections of soil water content
                        FgDepth: Heat flux plate depth, m
                        BulkDensity: Soil bulk density, kg/m3
                        OrganicContent: Soil organic content, fraction
                        SwsDefault
                        Constants for empirical corrections using log(sensor)
                            and exp(sensor) functions (SWC_a0, SWC_a1, SWC_b0,
                            SWC_b1, SWC_t, TDR_a0, TDR_a1, TDR_b0, TDR_b1,
                            TDR_t)
                        Variable and attributes lists (empSWCin, empSWCout,
                            empTDRin, empTDRout, linTDRin, SWCattr, TDRattr)
                [Output]:
                    Variable subset list for OzFlux output file
                [Variables]:
                    Variable names and parameters for:
                        Range check to set upper and lower rejection limits
                        Diurnal check to reject observations by time of day that
                            are outside specified standard deviation limits
                        Timestamps, slope, and offset for Linear correction
                [Plots]:
                    Variable lists for plot generation
            """
        logger.info(" Starting L3 processing ...")
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf) == 0:
            logger.info( " L3: no control file chosen")
            self.do_progress(text='Waiting for input ...')
            return
        infilename = pfp_io.get_infilenamefromcf(cf)
        if not pfp_utils.file_exists(infilename):
            self.do_progress(text='An error occurred, check the console ...')
            return
        self.ds2 = pfp_io.nc_read_series(infilename)
        if len(self.ds2.series.keys())==0:
            self.do_progress(text='An error occurred, check the console ...')
            del self.ds2
            return
        self.update_startenddate(str(self.ds2.series['DateTime']['Data'][0]),
                                 str(self.ds2.series['DateTime']['Data'][-1]))
        self.do_progress(text='Doing L3 QC & Corrections ...')
        self.ds3 = pfp_ls.l3qc(cf, self.ds2)
        # put a copy of the config file in self for later use
        self.cf = cf
        self.do_progress(text='Finished L3')
        txtstr = ' Finished L3: Standard processing for site: '
        txtstr = txtstr+self.ds3.globalattributes['site_name'].replace(' ','')
        logger.info(txtstr)
        self.do_progress(text='Saving L3 QC & Corrected NetCDF data ...')       # put up the progress message
        outfilename = pfp_io.get_outfilenamefromcf(cf)
        if len(outfilename)==0:
            self.do_progress(text='An error occurred, check the console ...')
            return
        ncFile = pfp_io.nc_open_write(outfilename)
        outputlist = pfp_io.get_outputlistfromcf(cf, 'nc')
        pfp_io.nc_write_series(ncFile, self.ds3, outputlist=outputlist)             # save the L3 data
        self.do_progress(text='Finished saving L3 QC & Corrected NetCDF data')  # tell the user we are done
        logger.info(' Finished saving L3 QC & Corrected NetCDF data')
        logger.info("")

    def do_l4qc(self):
        """
            Call pfp_ls.l4qc_gapfill function
            Performs L4 gap filling on L3 met data
            or
            Ingests L4 gap filled fluxes performed in external SOLO-ANN and c
                omputes daily sums
            Outputs L4 netCDF file to ncData folder
            Outputs L4 netCDF file to OzFlux folder

            ControlFiles:
                L4_year.txt
                or
                L4b.txt

            ControlFile contents (see ControlFile/Templates/L4.txt and
            ControlFile/Templates/L4b.txt for examples):
                [General]:
                    Python control parameters (SOLO)
                    Site characteristics parameters (Gap filling)
                [Files]:
                    L3 input file name and path (Gap filling)
                    L4 input file name and path (SOLO)
                    L4 output file name and ncData folder path (both)
                    L4 OzFlux output file name and OzFlux folder path
                [Variables]:
                    Variable subset list for OzFlux output file (where
                        available)
            """
        logger.info(" Starting L4 processing ...")
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        infilename = pfp_io.get_infilenamefromcf(cf)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        if not pfp_utils.file_exists(infilename): self.do_progress(text='An error occurred, check the console ...'); return
        ds3 = pfp_io.nc_read_series(infilename)
        if len(ds3.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del ds3; return
        ds3.globalattributes['controlfile_name'] = cf['controlfile_name']
        self.update_startenddate(str(ds3.series['DateTime']['Data'][0]),
                                 str(ds3.series['DateTime']['Data'][-1]))
        sitename = ds3.globalattributes['site_name']
        self.do_progress(text='Doing L4 gap filling drivers: '+sitename+' ...')
        if "Options" not in cf: cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        ds4 = pfp_ls.l4qc(cf,ds3)
        if ds4.returncodes["alternate"]=="quit" or ds4.returncodes["solo"]=="quit":
            self.do_progress(text='Quitting L4: '+sitename)
            logger.info(' Quitting L4: %s', sitename)
        else:
            self.do_progress(text='Finished L4: '+sitename)
            logger.info(' Finished L4: %s', sitename)
            self.do_progress(text='Saving L4 gap filled data ...')         # put up the progress message
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
            ncFile = pfp_io.nc_open_write(outfilename)
            outputlist = pfp_io.get_outputlistfromcf(cf,'nc')
            pfp_io.nc_write_series(ncFile,ds4,outputlist=outputlist)         # save the L4 data
            self.do_progress(text='Finished saving L4 gap filled data')    # tell the user we are done
            logger.info(' Finished saving L4 gap filled data')
        logger.info("")

    def do_l5qc(self):
        """
            Call pfp_ls.l5qc function to gap fill the fluxes.
        """
        logger.info("Starting L5 processing ...")
        cf = pfp_io.load_controlfile(path="controlfiles")
        if len(cf) == 0:
            self.do_progress(text="Waiting for input ...")
            return
        infilename = pfp_io.get_infilenamefromcf(cf)
        if len(infilename) == 0:
            self.do_progress(text="An error occurred, check the console ...")
            return
        if not pfp_utils.file_exists(infilename):
            self.do_progress(text="An error occurred, check the console ...")
            return
        ds4 = pfp_io.nc_read_series(infilename)
        if len(ds4.series.keys()) == 0:
            self.do_progress(text="An error occurred, check the console ...")
            del ds4
            return
        ds4.globalattributes["controlfile_name"] = cf["controlfile_name"]
        self.update_startenddate(str(ds4.series["DateTime"]["Data"][0]),
                                 str(ds4.series["DateTime"]["Data"][-1]))
        sitename = ds4.globalattributes["site_name"]
        self.do_progress(text="Doing L5 gap filling fluxes: "+sitename+" ...")
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        ds5 = pfp_ls.l5qc(cf, ds4)
        if ds5.returncodes["solo"] == "quit":
            self.do_progress(text="Quitting L5: "+sitename)
            logger.info(" Quitting L5: %s", sitename)
        else:
            self.do_progress(text="Finished L5: "+sitename)
            logger.info("Finished L5: %s", sitename)
            self.do_progress(text="Saving L5 gap filled data ...")           # put up the progress message
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            if len(outfilename)==0:
                self.do_progress(text="An error occurred, check the console ...")
                return
            ncFile = pfp_io.nc_open_write(outfilename)
            outputlist = pfp_io.get_outputlistfromcf(cf, "nc")
            pfp_io.nc_write_series(ncFile, ds5, outputlist=outputlist)         # save the L5 data
            self.do_progress(text="Finished saving L5 gap filled data")      # tell the user we are done
            logger.info("Finished saving L5 gap filled data")
        logger.info("")

    def do_l6qc(self):
        """
            Call pfp_ls.l6qc function to partition NEE into GPP and ER.
        """
        logger.info(" Starting L6 processing ...")
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        infilename = pfp_io.get_infilenamefromcf(cf)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        if not pfp_utils.file_exists(infilename): self.do_progress(text='An error occurred, check the console ...'); return
        ds5 = pfp_io.nc_read_series(infilename)
        if len(ds5.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del ds5; return
        ds5.globalattributes['controlfile_name'] = cf['controlfile_name']
        self.update_startenddate(str(ds5.series['DateTime']['Data'][0]),
                                 str(ds5.series['DateTime']['Data'][-1]))
        sitename = ds5.globalattributes['site_name']
        self.do_progress(text='Doing L6 partitioning: '+sitename+' ...')
        if "Options" not in cf: cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        ds6 = pfp_ls.l6qc(cf,ds5)
        self.do_progress(text='Finished L6: '+sitename)
        logger.info("Finished L6: %s", sitename)
        self.do_progress(text='Saving L6 partitioned data ...')           # put up the progress message
        outfilename = pfp_io.get_outfilenamefromcf(cf)
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        ncFile = pfp_io.nc_open_write(outfilename)
        outputlist = pfp_io.get_outputlistfromcf(cf,'nc')
        pfp_io.nc_write_series(ncFile,ds6,outputlist=outputlist)             # save the L6 data
        self.do_progress(text='Finished saving L6 partitioned data')      # tell the user we are done
        logger.info("Finished saving L6 partitioned data")
        logger.info("")

    def do_mpt(self,mode="standard"):
        """
        Calls pfp_mpt.mpt_main
        Calculate the u* threshold using the Moving Point Threshold (MPT) method.
        """
        logger.info(' Starting estimation u* threshold using MPT')
        self.do_progress(text='Estimating u* threshold using MPT ...')
        if mode=="standard":
            stdname = "controlfiles/standard/mpt.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                self.do_progress(text='Opening input file ...')
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if not os.path.exists(filename):
                    logger.info( " MPT: no input file chosen")
                    self.do_progress(text='Waiting for input ...')
                    return
                if "Files" not in dir(cf): cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                in_filename = ntpath.split(filename)[1]
                cf["Files"]["in_filename"] = in_filename
                cf["Files"]["out_filename"] = in_filename.replace(".nc","_MPT.xls")
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf)==0:
                    self.do_progress(text='Waiting for input ...')
                    return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        self.do_progress(text='Doing the u* threshold (MPT)')
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        pfp_mpt.mpt_main(cf)
        self.do_progress(text='Finished estimating u* threshold (MPT)')
        logger.info(' Finished estimating u* threshold (MPT)')
        logger.info("")

    def do_nc2ep_biomet(self):
        """ Calls pfp_io.ep_biomet_write_csv. """
        logger.info(' Starting conversion to EddyPro biomet file')
        self.do_progress(text='Load control file ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Converting nc to EddyPro biomet CSV ...')
        return_code = pfp_io.ep_biomet_write_csv(cf)
        if return_code==0:
            self.do_progress(text='An error occurred, check the console ...')
            return
        else:
            logger.info(' Finished conversion to EddyPro biomet format')
            self.do_progress(text='Finished conversion to EddyPro biomet format')
            logger.info("")

    def do_nc2csv_ecostress(self, mode="standard"):
        """ Calls pfp_io.write_csv_ecostress. """
        logger.info(' Starting conversion to ECOSTRESS file')
        if mode == "standard":
            stdname = "controlfiles/standard/nc2csv_ecostress.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if len(filename) == 0 or not os.path.exists(filename):
                    self.do_progress(text='Waiting for input ...')
                    return
                if "Files" not in dir(cf):
                    cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                cf["Files"]["in_filename"] = ntpath.split(filename)[1]
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf) == 0:
                    self.do_progress(text='Waiting for input ...')
                    return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf) == 0:
                self.do_progress(text='Waiting for input ...')
                return
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        return_code = pfp_io.write_csv_ecostress(cf)
        if return_code == 0:
            self.do_progress(text='An error occurred, check the console ...')
            return
        else:
            logger.info(' Finished conversion to EddyPro biomet format')
            self.do_progress(text='Finished conversion to EddyPro biomet format')
            logger.info("")

    def do_nc2fn(self):
        """ Calls pfp_io.fn_write_csv. """
        logger.info(' Starting conversion to FluxNet CSV file')
        self.do_progress(text='Load control file ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf) == 0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Converting nc to FluxNet CSV ...')
        pfp_io.fn_write_csv(cf)
        logger.info(' Finished conversion')
        self.do_progress(text='Finished conversion')
        logger.info("")

    def do_nc2reddyproc(self):
        """ Calls pfp_io.reddyproc_write_csv."""
        logger.info(' Starting conversion to REddyProc CSV file')
        self.do_progress(text="Load control file ...")
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Converting nc to REddyProc CSV ...')
        pfp_io.reddyproc_write_csv(cf)
        logger.info(' Finished conversion')
        self.do_progress(text='Finished conversion')
        logger.info("")

    def do_nc2smap(self):
        """ Calls pfp_io.smap_write_csv. """
        logger.info(' Starting conversion to SMAP CSV file')
        self.do_progress(text='Load control file ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf) == 0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Converting nc to SMAP CSV ...')
        pfp_io.smap_write_csv(cf)
        logger.info(' Finished conversion')
        self.do_progress(text='Finished conversion')
        logger.info("")

    def do_nc2xls(self):
        """ Calls pfp_io.nc_2xls. """
        logger.info(" Starting conversion to Excel file")
        self.do_progress(text="Choosing netCDF file ...")
        ncfilename = pfp_io.get_filename_dialog(path="../Sites", title="Choose a netCDF file")
        if len(ncfilename)==0:
            self.do_progress(text="Waiting for input ...")
            return
        self.do_progress(text="Converting netCDF file to Excel file")
        pfp_io.nc_2xls(ncfilename, outputlist=None)
        self.do_progress(text="Finished converting netCDF file")
        logger.info(" Finished converting netCDF file")
        logger.info("")

    def do_ncconcat(self):
        """
        Calls pfp_io.nc_concatenate
        """
        logger.info(' Starting concatenation of netCDF files')
        self.do_progress(text='Loading control file ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf) == 0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Concatenating files')
        pfp_io.nc_concatenate(cf)
        self.do_progress(text='Finished concatenating files')
        logger.info(' Finished concatenating files')
        logger.info("")

    def do_ncsplit(self):
        """
        Calls pfp_io.nc_split
        """
        logger.info(' Starting split of netCDF file')
        self.do_progress(text='Splitting file')
        pfp_io.nc_split()
        self.do_progress(text='Finished splitting file')
        logger.info(' Finished splitting file')
        logger.info("")

    def do_plotfingerprint(self,mode="standard"):
        """ Plot fingerprint"""
        logger.info(' Starting fingerprint plot')
        self.do_progress(text='Doing fingerprint plot ...')
        if mode=="standard":
            stdname = "controlfiles/standard/fingerprint.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if len(filename)==0 or not os.path.exists(filename):
                    self.do_progress(text='Waiting for input ...'); return
                if "Files" not in dir(cf): cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                cf["Files"]["in_filename"] = ntpath.split(filename)[1]
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        if "Options" not in cf:
            cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        self.do_progress(text='Plotting fingerprint ...')
        pfp_plot.plot_fingerprint(cf)
        self.do_progress(text='Finished plotting fingerprint')
        logger.info(' Finished plotting fingerprint')
        logger.info("")

    def do_plotfluxnet(self,mode="standard"):
        """ Plot FluxNet style time series of data."""
        self.do_progress(text='Doing FluxNet plots ...')
        if mode=="standard":
            stdname = "controlfiles/standard/fluxnet.txt"
            if os.path.exists(stdname):
                cf = pfp_io.get_controlfilecontents(stdname)
                filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose a netCDF file')
                if len(filename)==0 or not os.path.exists(filename):
                    self.do_progress(text='Waiting for input ...'); return
                if "Files" not in dir(cf): cf["Files"] = {}
                cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
                cf["Files"]["in_filename"] = ntpath.split(filename)[1]
            else:
                self.do_progress(text='Loading control file ...')
                cf = pfp_io.load_controlfile(path='controlfiles')
                if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        else:
            self.do_progress(text='Loading control file ...')
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Plotting FluxNet style plots ...')
        pfp_plot.plot_fluxnet(cf)
        self.do_progress(text='Finished FluxNet plotting')
        logger.info(' Finished FluxNet plotting')

    def do_plotquickcheck(self):
        """ Plot quickcheck"""
        self.do_progress(text='Loading control file ...')
        stdname = "controlfiles/standard/quickcheck.txt"
        if os.path.exists(stdname):
            cf = pfp_io.get_controlfilecontents(stdname)
            filename = pfp_io.get_filename_dialog(path='../Sites',title='Choose an input file')
            if len(filename)==0: self.do_progress(text='Waiting for input ...'); return
            if "Files" not in dir(cf): cf["Files"] = {}
            cf["Files"]["file_path"] = ntpath.split(filename)[0]+"/"
            cf["Files"]["in_filename"] = ntpath.split(filename)[1]
        else:
            cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Plotting quickcheck ...')
        pfp_plot.plot_quickcheck(cf)
        self.do_progress(text='Finished plotting quickcheck')
        logger.info(' Finished plotting quickcheck')

    def do_plotL1L2(self):
        """
            Plot L1 (raw) and L2 (QA/QC) data in blue and red, respectively

            Control File for do_l2qc function used.
            If L2 Control File not loaded, requires control file selection.
            """
        if 'ds1' not in dir(self) or 'ds2' not in dir(self):
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf) == 0:
                self.do_progress(text='Waiting for input ...')
                return
            l1filename = pfp_io.get_infilenamefromcf(cf)
            if not pfp_utils.file_exists(l1filename):
                self.do_progress(text='An error occurred, check the console ...')
                return
            ds1 = pfp_io.nc_read_series(l1filename)
            if len(ds1.series.keys()) == 0:
                self.do_progress(text='An error occurred, check the console ...')
                del ds1
                return
            l2filename = pfp_io.get_outfilenamefromcf(cf)
            ds2 = pfp_io.nc_read_series(l2filename)
            if len(ds2.series.keys()) == 0:
                self.do_progress(text='An error occurred, check the console ...')
                del ds2
                return
            self.update_startenddate(str(ds1.series['DateTime']['Data'][0]),
                                     str(ds1.series['DateTime']['Data'][-1]))
        else:
            ds1 = self.ds1
            ds2 = self.ds2
        self.do_progress(text='Plotting L1 & L2 QC ...')
        cfname = ds2.globalattributes['controlfile_name']
        cf = pfp_io.get_controlfilecontents(cfname)
        for nFig in cf['Plots'].keys():
            si = pfp_utils.GetDateIndex(ds1.series['DateTime']['Data'], self.plotstartEntry.get(),
                                      ts=ds1.globalattributes['time_step'], default=0, match='exact')
            ei = pfp_utils.GetDateIndex(ds1.series['DateTime']['Data'], self.plotendEntry.get(),
                                      ts=ds1.globalattributes['time_step'], default=-1, match='exact')
            plt_cf = cf['Plots'][str(nFig)]
            if 'Type' in plt_cf.keys():
                if str(plt_cf['Type']).lower() == 'xy':
                    self.do_progress(text='Plotting L1 and L2 XY ...')
                    pfp_plot.plotxy(cf, nFig, plt_cf, ds1, ds2, si, ei)
                else:
                    self.do_progress(text='Plotting L1 and L2 QC ...')
                    pfp_plot.plottimeseries(cf, nFig, ds1, ds2, si, ei)
            else:
                self.do_progress(text='Plotting L1 and L2 QC ...')
                pfp_plot.plottimeseries(cf, nFig, ds1, ds2, si, ei)
        self.do_progress(text='Finished plotting L1 and L2')
        logger.info(' Finished plotting L1 and L2, check the GUI')

    def do_plotL3L3(self):
        """
            Plot L3 (QA/QC and Corrected) data

            Control File for do_l3qc function used.
            If L3 Control File not loaded, requires control file selection.
            """
        if 'ds3' not in dir(self):
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf) == 0:
                self.do_progress(text='Waiting for input ...')
                return
            l3filename = pfp_io.get_outfilenamefromcf(cf)
            ds3 = pfp_io.nc_read_series(l3filename)
            if len(ds3.series.keys()) == 0:
                self.do_progress(text='An error occurred, check the console ...')
                del ds3
                return
            self.update_startenddate(str(ds3.series['DateTime']['Data'][0]),
                                     str(ds3.series['DateTime']['Data'][-1]))
        else:
            ds3 = self.ds3
        self.do_progress(text='Plotting L3 QC ...')
        cfname = ds3.globalattributes['controlfile_name']
        cf = pfp_io.get_controlfilecontents(cfname)
        for nFig in cf['Plots'].keys():
            si = pfp_utils.GetDateIndex(ds3.series['DateTime']['Data'], self.plotstartEntry.get(),
                                      ts=ds3.globalattributes['time_step'], default=0, match='exact')
            ei = pfp_utils.GetDateIndex(ds3.series['DateTime']['Data'], self.plotendEntry.get(),
                                      ts=ds3.globalattributes['time_step'], default=-1, match='exact')
            plt_cf = cf['Plots'][str(nFig)]
            if 'Type' in plt_cf.keys():
                if str(plt_cf['Type']).lower() =='xy':
                    self.do_progress(text='Plotting L3 XY ...')
                    pfp_plot.plotxy(cf, nFig, plt_cf, ds3, ds3, si, ei)
                else:
                    self.do_progress(text='Plotting L3 QC ...')
                    pfp_plot.plottimeseries(cf, nFig, ds3, ds3, si, ei)
            else:
                self.do_progress(text='Plotting L3 QC ...')
                pfp_plot.plottimeseries(cf, nFig, ds3, ds3, si, ei)
        self.do_progress(text='Finished plotting L3')
        logger.info(' Finished plotting L3, check the GUI')

    def do_plotL3L4(self):
        """
            Plot L3 (QA/QC and Corrected) and L4 (Gap Filled) data in blue and
                red, respectively

            Control File for do_l4qc function used.
            If L4 Control File not loaded, requires control file selection.
            """
        if 'ds3' not in locals() or 'ds4' not in locals():
            cf = pfp_io.load_controlfile(path='controlfiles')
            if len(cf) == 0:
                self.do_progress(text='Waiting for input ...')
                return
            l3filename = pfp_io.get_infilenamefromcf(cf)
            if not pfp_utils.file_exists(l3filename):
                self.do_progress(text='An error occurred, check the console ...')
                return
            ds3 = pfp_io.nc_read_series(l3filename)
            if len(ds3.series.keys())==0:
                self.do_progress(text='An error occurred, check the console ...')
                del ds3
                return
            l4filename = pfp_io.get_outfilenamefromcf(cf)
            ds4 = pfp_io.nc_read_series(l4filename)
            if len(ds4.series.keys()) == 0:
                self.do_progress(text='An error occurred, check the console ...')
                del ds4
                return
            self.update_startenddate(str(ds3.series['DateTime']['Data'][0]),
                                     str(ds3.series['DateTime']['Data'][-1]))
        self.do_progress(text='Plotting L3 and L4 QC ...')
        cfname = ds4.globalattributes['controlfile_name']
        cf = pfp_io.get_controlfilecontents(cfname)
        for nFig in cf['Plots'].keys():
            si = pfp_utils.GetDateIndex(ds3.series['DateTime']['Data'], self.plotstartEntry.get(),
                                      ts=ds3.globalattributes['time_step'], default=0, match='exact')
            ei = pfp_utils.GetDateIndex(ds3.series['DateTime']['Data'], self.plotendEntry.get(),
                                      ts=ds3.globalattributes['time_step'], default=-1, match='exact')
            pfp_plot.plottimeseries(cf, nFig, ds3, ds4, si, ei)
        self.do_progress(text='Finished plotting L4')
        logger.info(' Finished plotting L4, check the GUI')

    def do_plotL4L5(self):
        """
            Plot L4 (Gap filled) and L5 (Partitioned) data.
            """
        pass

    def do_plotL6_summary(self):
        """
            Plot L6 summary.
        """
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if "Options" not in cf: cf["Options"]={}
        cf["Options"]["call_mode"] = "interactive"
        l6filename = pfp_io.get_outfilenamefromcf(cf)
        if not pfp_utils.file_exists(l6filename): self.do_progress(text='An error occurred, check the console ...'); return
        ds6 = pfp_io.nc_read_series(l6filename)
        if len(ds6.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del ds6; return
        self.update_startenddate(str(ds6.series['DateTime']['Data'][0]),
                                 str(ds6.series['DateTime']['Data'][-1]))
        self.do_progress(text='Plotting L6 summary ...')
        pfp_gf.ImportSeries(cf,ds6)
        pfp_rp.L6_summary(cf,ds6)
        self.do_progress(text='Finished plotting L6 summary')
        logger.info(' Finished plotting L6 summary, check the GUI')

    def do_progress(self,text):
        """
            Update progress message in QC Data GUI
            """
        self.progress.destroy()
        self.progress = tk.Label(self.org_frame, text=text)
        self.progress.grid(row=8,column=0,columnspan=6,sticky="W")
        self.update()

    def do_quit(self):
        """
            Close plot windows and quit QC Data GUI
            """
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        logger.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        logger.info(' Quitting ...')
        self.quit()

    def do_savexL2(self):
        """
            Call nc2xl function
            Exports excel data from NetCDF file

            Outputs L2 Excel file containing Data and Flag worksheets
            """
        # check that we have a copy of the control file in self
        if not hasattr(self, "cf"):
            # load a control file if we don't
            self.cf = pfp_io.load_controlfile(path="controlfiles")
            if len(self.cf) == 0:
                logger.info( " No control file chosen")
                self.do_progress(text="Waiting for input ...")
                return
        self.do_progress(text="Exporting L2 NetCDF -> Excel ...")                     # put up the progress message
        # get the output filename
        outfilename = pfp_io.get_outfilenamefromcf(self.cf)
        # get the output list
        outputlist = pfp_io.get_outputlistfromcf(self.cf, "xl")
        pfp_io.nc_2xls(outfilename, outputlist=outputlist)
        # delete the control file attribute
        delattr(self, "cf")
        self.do_progress(text="Finished L2 Data Export")              # tell the user we are done
        logger.info(" Finished saving L2 data")

    def do_savexL3(self):
        """
            Call nc2xl function
            Exports excel data from NetCDF file

            Outputs L3 Excel file containing Data and Flag worksheets
            """
        # check that we have a copy of the control file in self
        if not hasattr(self, "cf"):
            # load a control file if we don't
            self.cf = pfp_io.load_controlfile(path="controlfiles")
            if len(self.cf) == 0:
                logger.info( " No control file chosen")
                self.do_progress(text="Waiting for input ...")
                return
        self.do_progress(text="Exporting L3 NetCDF -> Excel ...")                     # put up the progress message
        # get the output filename
        outfilename = pfp_io.get_outfilenamefromcf(self.cf)
        # get the output list
        outputlist = pfp_io.get_outputlistfromcf(self.cf, "xl")
        pfp_io.nc_2xls(outfilename, outputlist=outputlist)
        # delete the control file attribute
        delattr(self, "cf")
        self.do_progress(text="Finished L3 Data Export")              # tell the user we are done
        logger.info(" Finished saving L3 data")

    def do_xl2nc(self):
        """
        Calls pfp_io.xl2nc
        """
        logger.info(" Starting L1 processing ...")
        self.do_progress(text='Loading control file ...')
        cf = pfp_io.load_controlfile(path='controlfiles')
        if len(cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Reading Excel file & writing to netCDF')
        rcode = pfp_io.xl2nc(cf, "L1")
        if rcode==1:
            self.do_progress(text='Finished writing to netCDF ...')
            logger.info(' Finished writing to netCDF ...')
        else:
            self.do_progress(text='An error occurred, check the console ...')

    def update_startenddate(self, startstr, endstr):
        """
            Read start and end timestamps from data and report in QC Data GUI
            """
        self.filestartValue.destroy()
        self.fileendValue.destroy()
        self.filestartValue = tk.Label(self.org_frame, text=startstr)
        self.filestartValue.grid(row=3, column=0, columnspan=3)
        self.fileendValue = tk.Label(self.org_frame, text=endstr)
        self.fileendValue.grid(row=3, column=3, columnspan=3)
        self.update()

if __name__ == "__main__":
    qcGUI = qcgui(None)
    main_title = cfg.version_name+' Main GUI '+cfg.version_number
    qcGUI.title(main_title)
    qcGUI.mainloop()
    qcGUI.destroy()

    logger.info('PFP: All done')
