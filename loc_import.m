function [path_name, file_base, file_fmt] = loc_import
% Select .csv files to improt for location data from ezTrack
%   ZZ 04/10/19 modified from MIN1PIPE 

    [file_name_tmp, path_name] = uigetfile('*.*', 'Select first location tracking output file', 'MultiSelect', 'off');
    if ~iscell(file_name_tmp)
        file_name{1} = file_name_tmp;
    else
        file_name = file_name_tmp;
    end
    
    file_base = cell(1, length(file_name));
    file_fmt = cell(1, length(file_name));
    for i = 1: length(file_name)
        %%% first find out whether it is Inscopix or UCLA %%%
        if contains(file_name{i}, '.csv')
            file_fmt{i} = 'csv';
            ids = regexp(file_name{i}, '\d');
            ids = ids(end);
            file_base{i} = file_name{i}(1: ids - 1);
        elseif contains(file_name{i}, '.csvv')
            file_fmt{i} = 'csv';
            file_temp1 = find(file_name{i} == '-', 1, 'last');
            file_temp2 = find(file_name{i} == '.');
            if ~isempty(file_temp1)
                file_base{i} = file_name{i}(1: file_temp1 - 1);
            else
                file_base{i} = file_name{i}(1: file_temp2 - 1);
            end
        elseif contains(file_name{i}, '.tif')
            file_fmt{i} = 'tif';
            file_temp1 = find(file_name{i} == '-', 1, 'last');
            file_temp2 = find(file_name{i} == '.');
            if ~isempty(file_temp1)
                file_base{i} = file_name{i}(1: file_temp1 - 1);
            else
                file_base{i} = file_name{i}(1: file_temp2 - 1);
            end
        elseif contains(file_name{i}, '.mat')
            file_fmt{i} = 'mat';
            file_temp1 = find(file_name{i} == '-', 1, 'last');
            file_temp2 = find(file_name{i} == '.');
            if ~isempty(file_temp1)
                file_base{i} = file_name{i}(1: file_temp1 - 1);
            else
                file_base{i} = file_name{i}(1: file_temp2 - 1);
            end
        end
    end
end