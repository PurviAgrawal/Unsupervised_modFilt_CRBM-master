% This script passes the validation data spectrograms through CRBM (with
% its weight being fixed as the learned filters) and saves the hidden
% activation probabilities in function 'crbm_filtSelection_validation'
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

% Load the text file or list having paths of spectrograms of validation data
inList = 'validation_spec.txt';
[inSpec] = textread(inList, '%s');

total_feat = length(inSpec);
count_of_total_frames = 0;

feat = cell(1,total_feat);
for k = 1:total_feat
        
        load(inSpec{k,1});
                
        feat{1,k} = mel_spec_feat;  % for rate filtering of mel
%       feat{1,k} = mel_spec_feat';   % for scale filtering of mel
end

data1.y = feat;
clear feat

for k=1:total_feat
    data2.x(:,:,1,k) = data1.y(k);
end
clear data1

data.x = make_input_matrix_dev(data2.x); 
clear data2

% Compile mex files
make(0);

params = getparams;
params.verbose = 4;
for i = 1:size(data.x,4)
    data_one_utt = data.x(:,:,1,i);
    
    % For Rate filter: (else comment next line)
    data_one_utt = make_input_matrix_one_utt(data_one_utt);
    
    % For Scale filter: (else comment next line)
%     data_one_utt = make_input_matrix_one_utt(data_one_utt');
    
    [model ~] = crbm_filtSelection_validation(data_one_utt, params,i);
end
