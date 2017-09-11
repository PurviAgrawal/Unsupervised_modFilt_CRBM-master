function X=make_input_matrix_dev(input)
%
%   Written by: Purvi Agrawal
%   January, 2017
%   LEAP Lab, IISc

samples=floor(size(input,4));
% [H, W, C, N] = size(input);
C=size(input,3);
N=size(input,4);
% X = zeros(w(1),w(2),C,samples);
for i = 1:samples,
    im = randi(N);
    X{:,:,:,i} = input{:,:,1,im};
end
