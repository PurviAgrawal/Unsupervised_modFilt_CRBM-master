% This script loads the log mel spectrogram features of training data.
% For rate filter learning, the input to CRBM is sampled as temporal energy
% trajectories from spectrograms.
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

% Load the text file or list having paths of spectrograms of training data
inList = 'sourceFile_wav_spec.txt';
[inWav outSpec] = textread(inList, '%s %s');

%% For 2nd/3rd Rate filter learning (otherwise comment this section)

% % load previous learned 1st filter
% dataFile = 'data/model_rate/model_clean_mel_rate_15tap_1.mat';
% weights = load(dataFile);
% num_taps = length(weights.model.W);
% hh = abs(fft(weights.model.W,num_taps+1));
% scale_norm = 1/max(abs(hh));          %scaling factor
% hh = hh * scale_norm;
% % calculate residual
% hh=1-hh;
% F=[0:num_taps]/num_taps;
% w_rate1=firls(num_taps-1,F,hh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For 3rd Rate filter learning (otherwise comment this section)

% % load previous learned filters - 1st and 2nd
% dataFile = 'data/model_rate/model_clean_mel_rate_15tap_2.mat';
% weights = load(dataFile);
% num_taps = length(weights.model.W);
% hh = abs(fft(weights.model.W,num_taps+1));
% scale_norm = 1/max(abs(hh));          %scaling factor
% hh = hh * scale_norm;
% % calculate residual
% hh=1-hh;
% F=[0:num_taps]/num_taps;
% w_rate2=firls(num_taps-1,F,hh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
total_feat = length(outSpec);

feat = cell(1,total_feat);
for k = 1:total_feat
        
        load(outSpec{k,1});
        %%%% For 1st Rate filter learning (otherwise comment next line)
        feat{1,k} = mel_spec_feat'; % save in dimension of (no. of frames x no. of freq bands)
                
        %%%% For 2nd Rate filter learning   (otherwise comment next line)     
%         feat{1,k} = aud2cor_dataMod_rate(mel_spec_feat',w_rate1); % save in dimension of (no. of frames x no. of freq bands)
        
         %%%% For 3rd Rate filter learning (otherwise comment next 2 lines)
%         temp = aud2cor_dataMod_rate(mel_spec_feat',w_rate1);
%         feat{1,k} =  aud2cor_dataMod_rate(temp,w_rate2);  % save in dimension of (no. of frames x no. of freq bands)

end

data1.y = feat;
clear feat;

for k=1:total_feat
    data2.x(:,:,1,k) = data1.y(k); % to save the data in four dimensions
end

clear data1

size_of_min_frames=0;
for i=1:total_feat
    i
    frames = size(data2.x{:,:,:,i},1)
    if i==1
        size_of_min_frames = frames;
    end
    if frames <= size_of_min_frames
        size_of_min_frames = frames
    end
end
    

block_size = [size_of_min_frames 1];

% Sampling temporal energy trajectories from spectrograms for rate filter
% learning
data.x = make_input_matrix(data2.x, block_size);

clear data2

% Compile mex files
make(0);

params = getparams_rate;  % check the parameters of CRBM to be used
params.verbose = 4;
[model ~] = trainCRBM_Rate(data, params);
