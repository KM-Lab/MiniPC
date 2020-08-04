function [spkmap, gridprob, xedges, yedges, spkfn_fp, locs_dist_fp, spkmap_entry, gridcount] = binNFire_v3(px2cm, binsz, spkfn, locs_dist, spdreq, fps, frluX, fps_beh)
%% implement speed detection
% take speed and reject all activity while speed is < spdreq

spdminperframe = spdreq*px2cm/fps;  % use fps of miniscope, not fps of beh, because array pulled beh frames into ms frame times

% instanteous speed requirement applied at frame level
slow_frames = [];
j = 1;
for i = 1:length(locs_dist)
    if locs_dist(i,3) < spdminperframe
        slow_frames(j) = i;
        j = j+1;
    end
end

%% trim spk data to frames of interest and fast pass spk and location data
% trim spk data to frames of interset
spkfn_fp = [];
spkfn_fp = spkfn(:,frluX(:,1));

% remove both spike and location data from slow frames
spkfn_fp(:,slow_frames) = nan;
spkfn_fp = spkfn_fp(:,all(~isnan(spkfn_fp)));

locs_dist_fp = [];
locs_dist_fp = locs_dist;
locs_dist_fp(slow_frames,:) = nan;
locs_dist_fp(any(isnan(locs_dist_fp),2),:) = [];

% plot fast pass tracking
figure;
comet(locs_dist_fp(:,1), locs_dist_fp(:,2));
title('Location Tracking')


%% bin the space and get location probabilty map
binwidth = px2cm*binsz;
bnwidth2 = [binwidth binwidth];
xbnlm1 = min(locs_dist_fp(:,1))*0.75;
ybnlm1 = min(locs_dist_fp(:,2))*0.75;

xbnlm2 = max(locs_dist_fp(:,1))*1.1;
ybnlm2 = max(locs_dist_fp(:,2))*1.1;

% generate squre bins covering map and percent of time spent in bins
[gridprob, xedges, yedges] = histcounts2(locs_dist_fp(:,1), locs_dist_fp(:,2), 'xbinlimits', [0 xbnlm2], 'ybinlimits', [0 ybnlm2], 'BinWidth', bnwidth2, 'Normalization', 'probability');
[gridcount, ~, ~] = histcounts2(locs_dist_fp(:,1), locs_dist_fp(:,2), 'xbinlimits', [0 xbnlm2], 'ybinlimits', [0 ybnlm2], 'BinWidth', bnwidth2, 'Normalization', 'count');

gridprob(gridprob == 0) = nan;  % just for display purposes
figure;
heatmap(gridprob);
title('Location heatmap')
gridprob(isnan(gridprob)) = 0;


%% identify if a spk occurred; add it to spatial bin
nn = size(spkfn_fp,1);
xbin = 0;
ybin = 0;
this_x = [];
this_y = [];
last_x = 0;
last_y = 0;
spkmap = zeros(size(gridprob,1), size(gridprob,2), nn);
spkmap_entry = zeros(size(gridprob,1), size(gridprob,2));

for i = 1:nn
    for j = 1:size(spkfn_fp,2)
        if spkfn_fp(i,j) ~= 0
            x = locs_dist_fp(j,1);
            y = locs_dist_fp(j,2);
            for k = 1:size(gridprob,1)
                if xbin == 0
                     if (x >= xedges(k)) & (x < xedges(k+1))
                         xbin = k;
                     end
                end
            end
            for m = 1:(size(gridprob,2))
                if ybin == 0
                     if y >= yedges(m) & y < yedges(m+1)
                        ybin = m;
                     end
                end
            end 
            spkmap(xbin,ybin,i) = spkmap(xbin,ybin,i) + spkfn_fp(i,j);
            this_x = xbin;
            this_y = ybin;
            
            if i == 1
                if this_x ~= last_x | this_y ~= last_y 
                    spkmap_entry(xbin,ybin) = spkmap_entry(xbin,ybin)+1;
                end
            end
            last_x = this_x;
            last_y = this_y;
            xbin = 0;
            ybin = 0;
        end
    end
end

end