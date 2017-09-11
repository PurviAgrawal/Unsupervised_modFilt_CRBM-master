% This is an example to extract log mel spectrogram of audio files (which
% will be used for filter learning in next code)
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

addpath('melSpec_codes/')

% Load the text file or list containing paths of training data audio files and where the
% respective spectrograms to be saved, in 2 separate columns
inList = 'sourceFile_wav_spec.txt';
[inWav outSpec] = textread(inList, '%s %s');

for i = 1:length(inWav)
    [data fs]=audioread(inWav{i,1});

    % Bandpass filtering
    n = 7;          %order of the butterworth filter
    beginFreq = 200 / (fs/2);
    endFreq = 6500 / (fs/2);
    [b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
    data = filter(b, a, data);

    % LOG MEL SPECTROGRAM
    flen = 0.025*fs;       % frame length 25ms
    fhop = 0.010*fs;       % frame overlap 10ms

    fnum = floor((length(data) - flen) / fhop) + 1 ;
    req_len = (flen-1)*fhop + flen;          % Final number of required frames
    mel_spec_feat = melfcc(data,fs);
    mel_spec_feat = log(mel_spec_feat);
    mel_spec_feat = mel_spec_feat(:,1:fnum);

    save([outSpec{i,1}],'mel_spec_feat');       % save as feature vector

end

