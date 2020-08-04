function [SI, SI_rand, spkmap_rand, PC, smth_spkmap, smth_gridtime] = calcSI_v3_locshuf(spkfn_fp, locs_dist_fp, gridprob, nn, spkmap, numshuf, fps, xedges, yedges, fps_beh, gridcount)
%% what SI equation to use
SI_type = string('laure'); % options are 'laure' or 'david'

%% calculate spatial information (SI)
% lambda(x)*log
% spkgrd is spatial neural activity rate (spike count * bin occupancy)
% lambda_bar is overall mean firing rate
% lambda_i is firing rate for that bin
SI = zeros(nn,1);
gridprob_array = gridprob(:);
gridtime = gridcount./fps;
gridtime_array = gridtime(:);
sigma = 1; % in binsz units
    for i = 1:nn
        smth_spkmap(:,:,i) = imgaussfilt(spkmap(:,:,i),sigma);
    end
    smth_gridprob = imgaussfilt(gridprob,sigma);
    smth_gridtime = imgaussfilt(gridtime,sigma);
    smth_gridtime_array = smth_gridtime(:);
    smth_gridprob_array = smth_gridprob(:);
 
    if SI_type == 'david'
    % david method
        for i = 1:nn
           i_spkmap = smth_spkmap(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda(i) = sum(i_spkmap) / sum(smth_gridtime_array);
           lambda_x(i) = eps + (i_spkmap ./ (smth_gridtime_array + eps));
           p_x = smth_gridprob;
           I = sum(lambda_x(i) .* log2(lambda_x(i) ./ lambda(i)) .* p_x);
           SI(i) = I;
        end
    end
    
    if SI_type == 'laure'
    % laure method
        for i = 1:nn
           i_spkmap = smth_spkmap(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda_i = eps + (i_spkmap ./ (smth_gridtime_array + eps));
           lambda_bar(i) = eps + (sum(i_spkmap) / sum(smth_gridtime_array));
           p_x = smth_gridprob_array;
           I = sum((lambda_i ./ lambda_bar(i)) .* log2(lambda_i ./ lambda_bar(i)) .* p_x);
           SI(i) = I;
        end
    end


disp('Spatial information content calculated.');

%% next, shuffle location and bin spikes
disp('Now shuffling location, binning activity, and calculating SI');
SI_rand = zeros(nn, numshuf);
for perm = 1:numshuf
    if rem(perm,50) == 0
        disp(['Progress:  ' num2str(perm) ' out of ' num2str(numshuf) ' iterations complete.']);
    end
    % shuffle location
    rnd = randi([fps, length(locs_dist_fp)-fps], 1, 1); % trim off 1 second of location data from random sampling to minimize overlap
    locs_rand = [];
    locs_rand = locs_dist_fp(rnd:end,:);
    locs_rand = [locs_rand; locs_dist_fp(1:rnd-1,:)];
    
    % generate new spike map for each iteration
    xbin = 0;
    ybin = 0;
    spkmap_rand = zeros(size(gridprob,1), size(gridprob,2), nn);  
    for i = 1:nn
        for j = 1:size(spkfn_fp,2)
            if spkfn_fp(i,j) ~= 0
                x = locs_rand(j,1);
                y = locs_rand(j,2);
                for k = 1:size(gridprob,1)
                    if xbin == 0
                         if (x >= xedges(k)) && (x <= xedges(k+1))
                             xbin = k;
                         end
                    end
                end
                for m = 1:(size(gridprob,2))
                    if ybin == 0
                         if y >= yedges(m) && y <= yedges(m+1)
                            ybin = m;
                         end
                    end
                end
                %spkmap_rand(xbin,ybin,i) = spkmap_rand(xbin,ybin,i)+1;
                spkmap_rand(xbin,ybin,i) = spkmap_rand(xbin,ybin,i)+spkfn_fp(i,j);

                xbin = 0;
                ybin = 0;
            end
        end
    end
    
% calculate spatial information for this permutation    
    for i = 1:nn
        smth_spkmap_rand(:,:,i) = imgaussfilt(spkmap_rand(:,:,i),sigma);
    end

    if SI_type == 'david'
    % david method
        for i = 1:nn
           i_spkmap = smth_spkmap_rand(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda = sum(i_spkmap) / sum(smth_gridtime_array);
           lambda_x = eps + (i_spkmap ./ (smth_gridtime_array + eps));
           p_x = smth_gridprob_array;
           I = sum(lambda_x .* log2(lambda_x ./ lambda) .* p_x);
           SI_rand(i,perm) = I;
        end
    end


    if SI_type == 'laure'
    % laure method
        for i = 1:nn
           i_spkmap = smth_spkmap_rand(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda_i = eps + (i_spkmap ./ (smth_gridtime_array + eps));
           lambda_bar(i) = eps + (sum(i_spkmap) / sum(smth_gridtime_array));
           p_x = smth_gridprob_array;
           I = sum((lambda_i ./ lambda_bar(i)) .* log2(lambda_i ./ lambda_bar(i)) .* p_x);
           SI_rand(i,perm) = I;
        end
    end

end
disp('SI of randomized iterations complete');
%% determine if real SI reaches significance threshold (> 95 percentile)
PC = nan(nn,1);
for i = 1:nn
    cellofint = i;
    nless = sum(SI_rand(cellofint,:) < SI(cellofint));
    nequal = sum(SI_rand(cellofint,:) == SI(cellofint));
    centile = 100 * (nless + 0.5*nequal) / length(SI_rand(cellofint,:));
    if centile > 95
        PC(cellofint) = 1;
    else
        PC(cellofint) = 0;
    end
end

disp(['Number of cells with significant spatial information:  ' num2str(length((find(PC==1)))) ' out of ' num2str(length(PC))]);

end
