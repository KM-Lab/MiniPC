function [locs_dist, frluX, frlu, too_fast, locs] = ezlocTrack_v3_2TS(mscamnum, behavcamnum, dummybehavcamnum, frame_beh_sync, frame_ms_sync, frameStart, frameEnd, px2cm, fps_beh, fps_orig, fps)

%this version of ezlocTrack is for use when syncing a miniscope cam and a
%behavior cam recorded on different computers. Requires:  timestamp file
%for ms (with fake cam for syncing), timestamp file for real behavior cam,
%and manual prior identification of what frames on each are aligned

%% initialization
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');    % disables annoying warning message

[loc_path_name, loc_file_base, loc_file_fmt] = loc_import;
loc_path_name = char(loc_path_name);
loc_file_fmt = char(loc_file_fmt);
loc_file_base = char(loc_file_base);
loc_dirst = dir([loc_path_name, loc_file_base, '*', '.csv']);

[tfile, tpath] = uigetfile('.dat','Select miniscope timestamp file');  % file with timestamps of mscam and fake behcam
[t2file, t2path] = uigetfile('.dat','Select behavior timestamp file');  % file with timestamps of behcam

tbl2 = [];
locs = [];

%% concatenate location data

% convert location directory info to cell array to enable natural sorting
% of file names, uses natsort.m and notsortfiles.m (from matlab repo)
tmp_name = {loc_dirst.name};
tmp_name_sort = natsortfiles(tmp_name);
for i = 1:size(loc_dirst,1)
    filename = [loc_dirst(i).folder '\' char(tmp_name_sort(i))];
    tbl = [];
    tbl = readtable(filename);
    if i == 1
        tbl_all = tbl;
    else
        tbl_all = [tbl_all; tbl]; 
    end
end

locs = table2array(tbl_all(:,10:11));
disp('Finished concatenating location files');
%% align behavior cam and miniscope cam

tbl = readtable(strcat(tpath, tfile));
tbl = table2array(tbl);

tbl_beh = readtable(strcat(t2path, t2file));
tbl_beh = table2array(tbl_beh);

msframes = [];
dummybehavframes = [];
behavframes = [];

% separate frame table into miniscope table and dummy cam table
for i = 1:size(tbl,1)
    if tbl(i,1) == mscamnum
        msframes = [msframes; tbl(i,:)];
    elseif tbl(i,1) == dummybehavcamnum
        dummybehavframes = [dummybehavframes; tbl(i,:)];
    end
end
ds = fps_orig / fps; % temporal downsampling during processing
msframes = msframes(1:ds:end,:);    % adjusts for temporal ds

newnum = [1:size(msframes,1)];
newnum = newnum';
msframes(:,2) = newnum;
% get timestamp for frames that are synchronized and find offset between
% clocks
z = find(tbl(:,1) == dummybehavcamnum & tbl(:,2) == frame_ms_sync);
timesync_dummybeh = tbl(z,3);

z2 = find(tbl_beh(:,1) == behavcamnum & tbl_beh(:,2) == frame_beh_sync);
timesync_beh = tbl_beh(z2,3);

offset = timesync_beh - timesync_dummybeh;

% adjust behavior cam frames to align with sysclock of miniscope cam
behavframes = tbl_beh;
behavframes(:,3) = behavframes(:,3) - offset;


% for each miniscope frame, finds closest behavior frame 
match = [];
for i = 1:size(msframes,1)
    t = msframes(i,3);
    [~,ind] = min(abs(behavframes(:,3) - t));
     match(i) = ind;
end

disp('Aligned timepoints');
%% make a lookup table, col1 is miniscope frame, col2 is aligned behavior frame
% frlu stands for frame lookup
frlu = [];
frlu(:,1) = msframes(:,2);
frlu(:,2) = match;

% if no frame is specified for start of analysis, start at 1
if isempty(frameStart)
    frameStart = 1;
    % remove first frame (timestamp is junk)
    msframes = msframes([2:end],[1:3]);
    behavframes = behavframes([2:end],[1:3]);
end

% if no end frame is specified, go til end
if isempty(frameEnd)
    frameEnd = length(behavframes);
end

% only keep frames within start and end boundaries
[r, ~] = find(frlu(:,2) >= frameStart & frlu(:,2) <= frameEnd);
frluX = frlu(r,:);

% create array with locations of only the interested frames
locs_match = locs(frluX(:,2),1:2);

%% look for potential tracking errors
% note:  this feature isn't a  robust solution to poor tracking. Can only
% handle minor, quick jumps in tracking. Need to ensure tracking files are 
% of proper quality. 

spdlim = 100;  %in cm/s, set to a speed that clearly indicates tracking error
spdlim = spdlim*px2cm/fps;  %use fps of ms here, and not beh cam, because beh frames are pulled to match ms cam

    
dist = [];

% calculate distance moved each frame
for i = 2:length(locs_match)
    p1(1,:) = [locs_match(i,1) locs_match(i,2)];
    p1(2,:) = [locs_match(i-1,1) locs_match(i-1,2)];
    dist(i) = pdist(p1);
end

dist = dist';
too_fast = find(dist > spdlim);

if isempty(too_fast)
    display('No obvious tracking errors.');
else
    display('tracking error may have occurred. Removing offending timepoints');
end

% go through data and remove timepoints (and adjacent ones) for speeding
% frames
locs_dist = horzcat(locs_match, dist);
for i = 1:length(too_fast)
    locs_dist(too_fast(i)-1,:) = nan;
    locs_dist(too_fast(i),:) = nan;
    locs_dist(too_fast(i)+1,:) = nan;
    frluX(too_fast(i)-1,:) = nan;
    frluX(too_fast(i),:) = nan;
    frluX(too_fast(i)+1,:) = nan;
end

locs_dist(any(isnan(locs_dist),2),:) = [];
frluX(any(isnan(frluX),2),:) = [];
 
disp('Done analyzing tracking data.');

end