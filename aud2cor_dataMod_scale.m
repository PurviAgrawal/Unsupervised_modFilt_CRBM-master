function z = aud2cor_dataMod_scale(y,w)

% Filtering Log Mel Spectrogram
% w is a set of weights (flipped version of filter) obtained from CRBM
% here w is scale filter

% Sriram Ganapathy
% Leap Labs, IISc 
% 03-09-2016

[N,M] = size(y) ;

meany = mean(y(:));
% tempMean = repmat(meany,N,1);

y1 = y - meany ; % remove data mean temporally

w = w(:);

f = flipud(w) ; % make it as a filter

 F = repmat(abs(fft(f,N)),1,M); % Filter absolute FFT (The phase response has delaying effect)
 Y = fft(y1,N,1);   % Apply fft on temporal dimension
 Z = F.*Y ;
 z = ifft (Z,N,1);
 
 z=z';
 

