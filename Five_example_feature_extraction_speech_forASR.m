% This is an example to extract proposed features of audio files using
% the learned and selected filters through CRBM
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
clear all; clc

addpath('mel_spectrogram/')

% Load the text file or list containing paths of training data audio files and where the
% respective filtered spectrograms to be saved, in 2 separate columns
inList = 'sourceFile_wav_finalFeat.txt';
[inWav outSpec] = textread(inList, '%s %s');

for i = 1:length(inWav)
    [data fs]=audioread(inWav{i,1});
    % Bandpass filtering
    n = 7; %order of the butterworth filter
    beginFreq = 200 / (fs/2);
    endFreq = 6500 / (fs/2);
    [b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
    data = filter(b, a, data);

    % LOG MEL SPECTROGRAM
    flen=0.025*fs;  % frame length 25ms
    fhop=0.010*fs;  % frame overlap 10ms

    fnum = floor((length(data) - flen) / fhop) + 1 ;
    req_len = (flen-1)*fhop + flen;  % Final number of required frames
    mel_spec_feat = melfcc(data,fs);
    mel_spec_feat = log(mel_spec_feat);
    mel_spec_feat = mel_spec_feat(:,1:fnum);    % Log mel spectrogram (no. of bands x no. of frames) which is to be filtered

    %% RATE FILTERS
    dataFile1 = 'data/model_rate/model_clean_mel_rate_15tap_1.mat';
    m1 =  load (dataFile1) ;
    w1 = m1.model.W ; w1 = w1(:);

    dataFile2 = 'data/model_rate/model_clean_mel_rate_15tap_2.mat';
    m2 =  load (dataFile2) ;
    w2 = m2.model.W ; w2 = w2(:);

    % weights_rate = [w1 w2];
    weights_rate = [w2];      % select filter

    %% SCALE FILTERS
    dataFile3 = 'data/model_scale/model_clean_mel_scale_8tap_1.mat';
    m3 =  load (dataFile3) ;
    w3 = m3.model.W ; w3 = w3(:);

    dataFile4 = 'data/model_scale/model_clean_mel_scale_8tap_2.mat';
    m4 =  load (dataFile4) ;
    w4 = m4.model.W ; w4 = w4(:);

    %weights_scale = [w3];
    weights_scale=[w3 w4];   % select filter

    filt_Mel_Spec = cell (2,1);  % ACC. TO NO. OF FILTERS *****
    feat = [];

    %% RATE AND SCALE FILTERING- DATA DRIVEN

    for kk = 1 : length(filt_Mel_Spec)
        temp = aud2cor_dataMod_rate(mel_spec_feat',weights_rate(:,1));
        filt_Mel_Spec{kk} = aud2cor_dataMod_scale(temp',weights_scale(:,kk));
        feat = [ feat filt_Mel_Spec{kk}]; % concatenate the two filtered spectrograms
    end

    feat = feat'; % feat is the proposed feature for each input speech file

    save([outSpec{i,1}],'feat');       % save as feature vector

end

