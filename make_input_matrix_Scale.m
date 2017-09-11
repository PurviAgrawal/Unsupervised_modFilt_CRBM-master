function X=make_input_matrix_Scale(input, w, total_frames_all_data)
%
%   Written by: Purvi Agrawal
%   January, 2017
%   LEAP Lab, IISc

% samples = 100000;
% samples = 750000;
samples=floor(total_frames_all_data/2);
% [H, W, C, N] = size(input);
C=size(input,3);
N=size(input,4);
% patch_size=[113,21];
X = zeros(w(1),w(2),C,samples);
for i = 1:samples,
    im = randi(N);
    [H W]=size(input{:,:,1,im});
    temp_data=input{:,:,1,im};
    x = randi(W-w(1)+1)-1; y = randi(H-w(2)+1)-1;
    patch = temp_data(y+1:y+w(2), x+1:x+w(1));
    X(:,:,1,i) = patch;
    clear temp_data
    clear patch
end