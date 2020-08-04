# MiniPC
Place cell identification and analysis from miniscope recordings

Example dataset can be found [here](https://drive.google.com/drive/folders/1EBUdYbQ764NmUG4B9_i2_M5lvZpm_5Rh?usp=sharing).

The goal of Mini PC is to take processed miniscope data (from MIN1PIPE pipeline; [Lu et al., 2018](https://github.com/JinghaoLu/MIN1PIPE)) and location data (from ezTrack; [Pennington et al., 2019](https://github.com/DeniseCaiLab/ezTrack)) and identify place cells based on spatial information.

Run from script MiniPC.m

This software was written in part to solve the problem of aligning the timestamps for the miniscope recording with a seperate timestamp for a behavior video recording (eg when the behavior video has to be recorded on a seperate computer).

To synchronize the two timestamps, a synchronization cue (eg a light flash) has to be manually identified in both miniscope-behavior camera (called 'dummybehavcam' in the sofrware) with the actual behavior camera (called 'behavcam'). These synchronizing frames are to be entered into frame_beh_sync for the behavcam, and frame_ms_sync for the dummybehavcam. This allows an offset between the two timestamps to be identified, thus putting the two timestamps in the same relative timing.


**********
INPUTS

Manually entered in settings of MiniPC.m:  
- cam numbers
- frame synchronization numbers
- size calibrations
- behavior camera fps

Prompted for file input:  
- processed calcium file from MIN1PIPE
- timestamp file for msCam (from miniscope DAQ recording) 
- timestamp file for behaviorCam 
- animal location data (from ezTrack output .csv files assuming two object locations, required for formatting)

OUTPUTS
- Folder containing activity .fig rate maps ("spiking maps") for all place cells
- .mat file saved with various outputs and params

Key output variables:

id:  ID of mouse,

PC:  binary of whether or not cell has significant spatial information (ie a place cell),

SI:  spatial information (in bits) for each cell,

SI_rand: SI calculations for shuffled data,

spkmap:  neural activity rate map (sum of activity in each bin, normalized for time in bin),

smth_spkamp:  spkmap with gaussian filter applied,

locs_dist: matrix with columns for X coordinate, Y coordinate, and distance travelled for each frame,

locs_dist_fp (fp = fast pass):  locs_dist matrix with only portions of interest included, also speed filtered to remove frames in which animal was moving slower than spdreq,

spkfn_fp: spkfn data from processed calcium file, also fast passed,

gridprob:  heatmap of space with percent time spent in that bin,

frluX (frame lookup edit):  lookup table with miniscope frame number on the left, corresponding behavior frame on the right, edited to only include frames of interest (fp),

PC_info:  compiled information on place cells; row 1) cell ID, row 2) place field size, row 3) place field angle (in polar coordinates),

per_binoc:  percent of bins occupied in the arena (*specific for round arena)
