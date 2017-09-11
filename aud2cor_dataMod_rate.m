function z = aud2cor_dataMod_rate(y,w)

% Filtering Mel Spectrogram
% w is a set of weights (flipped version of filter) obtained from RBM
% here w is rate filter

% Sriram Ganapathy
% Leap Labs, IISc 
% 03-09-2016

[N,M] = size(y) ;

meany = mean(y,1);
tempMean = repmat(meany,N,1);

y1 = y - tempMean ; % remove data mean temporally
y1 = y;
w = w(:);

f = flipud(w) ; % make it as a filter

F = repmat(abs(fft(f,N)),1,M); % Filter absolute FFT (The phase response has delaying effect)
Y = fft(y1,N,1);   % Apply fft on temporal dimension
Z = F.*Y ;
z = ifft (Z,N,1);

