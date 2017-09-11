function X = make_input_matrix_one_utt(input)
%
%   Written by: Purvi Agrawal
%   January, 2017
%   LEAP Lab, IISc

samples = size(input{1,1}(:,1),1);  % no. of bands
w = [size(input{1,1}(:,:),2) 1];  % size of each band (no. of frames) x 1
X = zeros(w(1),w(2),1,samples);  % size of each band (no. of frames) x 1 x 1 x no. of bands
for i = 1:samples
    temp = input{1,1}(i,:);
    X(:,:,1,i) = temp(:);
end
