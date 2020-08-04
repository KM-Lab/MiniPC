% MiniPC (MINIscope Place Cells)   
% updated 2020.08.03

% for aligning animal location and calcium imaging data, then identifying
% place cells and visualzing place fields

% PREREQUISITES
% msCam_data_processed.mat output from MIN1PIPE (Lu et al., 2018, Cell Rep)
% behavior tracking done with ezTrack (Pennington et al., 2019, Sci Rep)
% using two objects (req'd for data formatting)
% using circular arena for behavior (if not, need to edit hardcode in
% analyzePCs function)

% INPUTS
% manual in settings:  cam numbers, frame synchronization numbers, size
% calibrations, behavior camera fps
% prompted:  processed calcium file from MIN1PIPE, timestamp file for msCam 
% and animalCam, animal location data from ezTrack output

% OUTPUTS
% Folder in MiniPC directory containing place fields and a .mat file with a 
% list of place cells, calcium data, location data, speed filtered data,
% spatial information, params, and other calculcations

%% settings
clear;
mscamnum = 0;   % camera number recording miniscope
dummybehavcamnum = 1;  % camera number on timestamp file from miniscope computer; used to synchronize other camera
behavcamnum = 0;    % camera number recording behavior

dorawCa2 = true; % choose whether to use deconvolved or raw calcium transients

frame_beh_sync = 662;  % synced frame on behavior cam
frame_ms_sync = 326;    % synced frame on ms_beh cam (dummybehavcam) 

frameStart = 1050; % behavior frame to begin  ([] = beginning)
frameEnd = 22000;  % behavior frame to end ([] = end)

px2cm = 8; % how many pixels in a centimeter
arenasz = 50; % arena size (assumed circular) in cm
binsz = 2.5;  % how many cm x cm in each spatial bin (binsz x binsz cm^2 bins)
spdreq = 2.5;   % min speed (in cm/s) for data to be included
fps_beh = 30; %fps of the behavior camera
numshuf = 500;  % number of shuffles to do for shuffled analysis

ID = '80v';  % name of mouse
id = strcat(ID,'_test'); % includes session info

dosave = true;  % true or false

%% load calcium data and frame timestamps
[cafile, capath] = uigetfile('Select processed calcium data');
if dorawCa2 == 1
    load(strcat(capath, cafile), 'sigfn');
    spkfn = sigfn; % code raw calcium data as spk data for use in analysis
    load(strcat(capath, cafile), 'Params');
    fps_orig = Params.Fsi;
    fps = Params.Fsi_new;
    nn = size(sigfn,1);
else
    load(strcat(capath, cafile), 'spkfn');
    load(strcat(capath, cafile), 'sigfn');
    load(strcat(capath, cafile), 'Params');
    fps_orig = Params.Fsi;
    fps = Params.Fsi_new;
    nn = size(spkfn,1);
end
%% location tracking concatenation (and minor error correction)
    [locs_dist, frluX, frlu, too_fast, locs] = ezlocTrack_v3_2TS(mscamnum, behavcamnum, dummybehavcamnum, frame_beh_sync, frame_ms_sync, frameStart, frameEnd, px2cm, fps_beh, fps_orig, fps);

    
%% bin pixel space into cm and make firing maps
    [spkmap, gridprob, xedges, yedges, spkfn_fp, locs_dist_fp, spkmap_entry, gridcount] = binNFire_v3(px2cm, binsz, spkfn, locs_dist, spdreq, fps, frluX, fps_beh);

%% calculate spatial information
    [SI, SI_rand, spkmap_rand, PC, smth_spkmap, smth_gridtime] = calcSI_v3_locshuf(spkfn_fp, locs_dist_fp, gridprob, nn, spkmap, numshuf, fps, xedges, yedges, fps_beh, gridcount);

%% analyze PCs
[PC_info, pc_ratemap, per_binoc] = analyzePCs (PC, spkmap, id, fps, gridprob, spkfn_fp, binsz, arenasz);

%% save
if dosave
    filename = strcat(['behPC_results_' id '.mat']);
    Params = struct;
    Params.mscamnum = mscamnum;
    Params.dummybehavcamnum = dummybehavcamnum;
    Params.behavcamnum = behavcamnum;
    Params.dorawCa2 = dorawCa2;
    Params.frame_beh_sync = frame_beh_sync;
    Params.frame_ms_sync = frame_ms_sync;
    Params.frameStart = frameStart;
    Params.frameEnd = frameEnd;
    Params.px2cm = px2cm;
    Params.binsz = binsz;
    Params.spdreq = spdreq;
    Params.fps_beh = fps_beh;
    Params.numshuf = numshuf;
    
    save([pwd '\' id '\' filename], 'Params', 'id', 'PC', 'SI', 'SI_rand', 'spkmap', 'smth_spkmap', 'spkfn_fp', 'spkfn', 'sigfn', 'gridprob', 'locs_dist_fp', 'locs_dist', 'frluX', 'PC_info', 'pc_ratemap', 'fps', 'per_binoc', 'spkmap_entry');
    disp('Analysis complete.');
end
