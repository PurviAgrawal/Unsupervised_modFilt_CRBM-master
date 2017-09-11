% This script computes the average of hidden activation probability values for 
% the validation data spectrograms (this value is to be compared for
% different rate and scale filters)
%
%
% Copyright (C) 2017 Purvi Agrawal & Sriram Ganapathy
%
%   January, 2017
%   LEAP Lab, IISc
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

source_path = 'data/model_hid_prob/clean_mel_15tap_rate1'; % Load the filter directory whose average hidden activation is to be computed (1 out of 6)
data_names=dir(source_path);
data_folder=char(data_names.name);
no_of_data_files=size(data_folder,1);

error_utt = zeros(1,no_of_data_files-2);
hid_prob_all = cell(1,1,no_of_data_files-2);
hid_prob_mean_each_utt = zeros(1,no_of_data_files-2);

% FOR SCALE FILTERING ONLY  (optional)
%hid_prob_mean_each_band = zeros(33,no_of_data_files-2); % for mel


for i = 3: no_of_data_files
    i
    file = strcat(source_path,'/',data_folder(i,:));
    load(file);
    error_utt(i-2) = errsum;
    temp = zeros(size(poshidprobs,1),size(poshidprobs,4));
    for k = 1: size(poshidprobs,4)
        temp(:,k) = poshidprobs(:,1,1,k);
    end
    hid_prob_all {1,1,i-2} = temp;
    hid_prob_mean_each_utt(i-2) = mean(mean(temp,1),2);
   % hid_prob_mean_each_band(:,i-2) = mean(temp,2); % optional
end

err_avg = mean(error_utt);
hid_prob_avg_all_utt = mean(hid_prob_mean_each_utt);  % the final avg. value which is to be noted and compared 

