function [SI, SI_rand, spkmap_rand, PC, smth_spkmap, smth_gridtime] = calcSI_v3(spkfn_fp, locs_dist_fp, gridprob, nn, spkmap, numshuf, fps, xedges, yedges, fps_beh, gridcount, smooth)
%% what SI equation to use
SI_type = string('rr'); % options are 'rr' (rondi-reig) or 'redish'

%% calculate spatial information (SI)
% lambda(x)*log
% spkgrd is spatial neural activity rate (spike count * bin occupancy)
% lambda_bar is overall mean firing rate
% lambda_i is firing rate for that bin
SI = zeros(nn,1);
gridprob_array = gridprob(:);
gridtime = gridcount./fps;
gridtime_array = gridtime(:);
smth_spkmap = [];
smth_gridprob = [];
if smooth == true;
    
    for i = 1:nn
        smth_spkmap(:,:,i) = imgaussfilt(spkmap(:,:,i),2);
    end
    smth_gridprob = imgaussfilt(gridprob,2);
    smth_gridtime = imgaussfilt(gridtime,2);
    smth_gridtime_array = smth_gridtime(:);
    smth_gridprob_array = smth_gridprob(:);
 
    if SI_type == string('redish')
    % redish method
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


    if SI_type == string('rr')
    % rochefort method
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


else
    if SI_type == string('redish')
    % redish method
        for i = 1:nn
            i_spkmap = spkmap(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda = sum(i_spkmap) / sum(gridtime_array);
           lambda_x = eps + (i_spkmap ./ (gridtime_array + eps));
           p_x = gridprob;
           I = sum(lambda_x .* log2(lambda_x ./ lambda) .* p_x);
           SI(i) = I;
        end
    end


    if SI_type == string('rr')
    % rochefort method
        for i = 1:nn
           i_spkmap = spkmap(:,:,i);
           i_spkmap = i_spkmap(:);

           lambda_i = eps + (i_spkmap ./ (gridtime_array + eps));
           lambda_bar(i) = eps + (sum(i_spkmap) / sum(gridtime_array));
            p_x = gridprob_array;
           I = sum((lambda_i ./ lambda_bar(i)) .* log2(lambda_i ./ lambda_bar(i)) .* p_x);
           SI(i) = I;
        end
    end

end
disp('Spatial information content calculated.');

%% shuffle trials  %%
disp('Now shuffling data.');

% Shuffles timing of spike (via permutation of frame #s)

spk_shuf = zeros(size(spkfn_fp,1), size(spkfn_fp,2), numshuf);
for shuf = 1:numshuf
    if rem(shuf,50) == 0
        disp(['Shuffle progress:  ' num2str(shuf) ' out of ' num2str(numshuf) ' complete.']);
    end
    spkfn_rand = zeros(size(spkfn_fp));
    % shuffle spike timing
    for i = 1:nn
        spks = spkfn_fp(i,:);
        spk_rand = spks(randperm(length(spks)));
        spkfn_rand(i,:) = spk_rand;
    end
    spk_shuf(:,:,shuf) = spkfn_rand;
end

disp('Spike timing shuffled');


%% next, bin shuffled spikes and calculate SI. Iterate over shuffles.
disp('Now calculating spatial information of shuffled data');
SI_rand = zeros(nn, numshuf);
for perm = 1:numshuf
    if rem(perm,50) == 0
        disp(['Progress:  ' num2str(perm) ' out of ' num2str(numshuf) ' iterations complete.']);
    end
    % assign current permutation of data to spk handle
    spk_shuf_i = spk_shuf(:,:,perm);
    % generate new spike map for each iteration
    nn = size(spk_shuf_i,1);
    xbin = 0;
    ybin = 0;
    spkmap_rand = zeros(size(gridprob,1), size(gridprob,2), nn);  
    for i = 1:nn
        for j = 1:size(spk_shuf_i,2)
            if spk_shuf_i(i,j) ~= 0
                x = locs_dist_fp(j,1);
                y = locs_dist_fp(j,2);
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
                spkmap_rand(xbin,ybin,i) = spkmap_rand(xbin,ybin,i)+1;
                xbin = 0;
                ybin = 0;
            end
        end
    end
    
  % calculate spatial information for this permutation
    if smooth == true;
    
        for i = 1:nn
            smth_spkmap_rand(:,:,i) = imgaussfilt(spkmap_rand(:,:,i),2);
        end

        if SI_type == string('redish')
        % redish method
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


        if SI_type == string('rr')
        % rochefort method
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


    else
        if SI_type == string('rr')  
          % rochefort method
            for i = 1:nn
              i_spkmap = spkmap_rand(:,:,i);
              i_spkmap = i_spkmap(:);

              lambda_i = eps + (i_spkmap ./ (gridtime_array + eps));
              lambda = eps + (sum(i_spkmap) / sum(gridtime_array));
               p_x = gridprob_array;
              I = sum((lambda_i ./ lambda) .* log2(lambda_i ./ lambda) .* gridtime);
             SI_rand(i,perm) = I;
            end
        end

        if SI_type == string('redish')
        % redish method
         for i = 1:nn
               i_spkmap = spkmap_rand(:,:,i);
               i_spkmap = i_spkmap(:);

               lambda = sum(i_spkmap) / sum(gridtime_array);
               lambda_x = eps + (i_spkmap ./ (gridtime_array +eps));
               p_x = gridprob_array;
               I = sum(lambda_x .* log2(lambda_x ./ lambda) .* p_x);
               SI_rand(i,perm) = I;
        end    
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
    if centile >= 95
        PC(cellofint) = 1;
    else
        PC(cellofint) = 0;
    end
end

disp(['Number of cells with significant spatial information:  ' num2str(length((find(PC==1)))) ' out of ' num2str(length(PC))]);

end
%%
% spk_n = spkmap(:,:,i);
% figure;
% heatmap(spk_n);

% 
% figure;
% for i = 1:nn
%     spk_n = spkmap(:,:,i);
%     heatmap(spk_n);
%     pause(1);
%     clf;
% end