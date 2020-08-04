function [PC_info, pc_ratemap, per_binoc] = analyzePCs (PC, spkmap, id, fps, gridprob, spkfn_fp, binsz, arenasz)

% make dirs to save figures and eventually variables
mkdir([id]);
mkdir([pwd '\' id], 'PlaceFields');
pcs = find(PC==1);
pc_bintime = gridprob * length(spkfn_fp) / fps;
pc_bintime_smth = imgaussfilt(pc_bintime,1);
pc_bintime_smth(gridprob==0) = 0;
    
pc_spkmap = zeros(size(spkmap,1),size(spkmap,2),length(pcs));
pc_ratemap = zeros(size(spkmap,1),size(spkmap,2),length(pcs));

for i = 1:length(pcs)
    pc_spkmap(:,:,i) = spkmap(:,:,pcs(i));
    pc_spkmap(:,:,i) = imgaussfilt(pc_spkmap(:,:,i),1);
end

% create rate maps for each PC (spks normalized to time)
for i = 1:length(pcs)
    i_pc_ratemap = pc_spkmap(:,:,i) ./ pc_bintime_smth;
    i_pc_ratemap(gridprob==0) = nan;
    pc_ratemap(:,:,i) = i_pc_ratemap;
end

% create and save rate map figs for all PCs
for i = 1:length(pcs)
    this_name = ['Cell' num2str(pcs(i)) 'RateMap.fig'];
    figure;
    h = heatmap(pc_ratemap(:,:,i));
    h.Colormap = jet;
    title([id, this_name]);
    savefig([pwd '\' id '\' 'PlaceFields\' this_name]);
end


% find size of place field and angle of PF
for i = 1:length(pcs)    
    % find size of PF (connected bins > 50% peak firing rate)
    i_ratemap = pc_ratemap(:,:,i);
    pk_rate = max(max(i_ratemap));
    i_ratemap(i_ratemap <= (pk_rate*0.5)) = 0;
    i_ratemap(isnan(i_ratemap)) = 0;
    [L, ~] = bwlabel(i_ratemap,8);
 
    % find the labeling associated with the peak rate
    [r, c] = find(i_ratemap == pk_rate);
    n = L(r,c);
    
    pf_size = length(find(L==n));

    % find angle of PF (area of peak firing rate in polar coordinates)
    new_coords_x = c-(size(i_ratemap,1)/2);
    new_coords_y = (size(i_ratemap,2)/2) - r;
    [pf_angle, ~] = cart2pol(new_coords_x, new_coords_y);
    
    % PC_info matrix:  row1: cellID; row2: PF size; row3: PF angle
    PC_info(1,i) = pcs(i);  
    PC_info(2,i) = pf_size;
    PC_info(3,i) = pf_angle;

end

% determine what percent of spatial bins were visited

% HARD CODE FOR  SHAPE (binarea uses circle formula)
% NEEDS TO BE CHANGED IF USING DIFFERENT ARENA
arenasz = arenasz / binsz;
binarea = round(pi*(arenasz/2)^2);

gp = gridprob;
gp(gridprob ~= 0) = 1;
oc_bins = find(gp(:,:)==1);
per_binoc = length(oc_bins) / binarea;

disp([int2str(per_binoc*100) '% of bins visited']);
close all

end