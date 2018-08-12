version_name = "PyFluxPro"
version_number = "V0.1.5"
# V0.1.5 - change source file names from qc*.py to pfp_*.py
# V0.1.4 - implement MPT and MDS
#          - implementation of u* threshold detection using the
#            Moving Point Threshold (MPT) technique via the
#            FluxNet C code
#          - implementation of Marginal Distribution Sampling (MDS)
#            gap filling via the FluxNet C code
#            - note that minor change made to common.c around
#              line 535 to fix bug that caused Calperum to fail
# V0.1.3 - "make it faster" version
#        - at some stage, about the time OFQC became PFP, PFP
#          became very slow, especially monthly summaries at
#          L6.  The suspect was the use of pfp_utils.GetVariables()
#          which I've been using to replace GetSeriesasMA().
#          Some quick profiling in Jupyter showed that
#          GetSeriesasMA() took ~750 us to produce the data
#          compared to ~350 ms for GetVariables(), that's a
#          factor of 500 slower!
#        - the culprit turned out to be the use of numpy.array()
#          to convert the datetime (originally a list) to an
#          array (~85 ms) plus some unnecessary use of
#          pfp_utils.get_start_index() and get_end_index().
#        - the fix was to re-write PFP so that datetime is stored
#          in the data structure as an array rather than a list plus
#          some changes to GetDateIndex(), FindIndicesOfBInA()
#          and sundry other routines.
# V0.1.2 - version produced during resolution of differences
#          between OFQC V2.9.6 and PFP V0.1.1 found by Craig
#          McFarlane at Great Western Woodlands:
#          - auto-complete at L4 was effectively disabled due
#            to a line inserted during debugging being left in
#            place in pfp_gf.gfalternate_get_correcteddata().
#          - reinstated original line to re-enable auto-complete
# V0.1.1 - cumulative changes up to 13/07/2017
#        - major change to pfp_utils.CreateSeries()
#          - deprecated "FList=" argument
#          - made all arguments compulsory
#        - QC flag generated for each series immediately
#          prior to pfp_utils.CreateSeries() based solely
#          on series data mask
# V0.1.0 - copy of OzFluxQC V2.9.6f
