function params = getparams(method)
% GETPARAMS  Get default params for trainCRBM
%
%   See also TRAINCRBM
%
%   Written by: Peng Qi, Sep 27, 2012
%
%   Original code can be found at : http://qipeng.me/software/convolutional-rbm.html
%   https://github.com/qipeng/convolutionalRBM.m. For research / personal
%   purposes only.

%% Model parameters
params.nmap = 1;
params.szFilter = [15 1];  % [15 1] for rate, [8 1] for scale filter learning
params.szPool = 1;
params.method = 'CD';

if (nargin > 0)
    if strcmp(method, 'PCD'),
        params.method = 'PCD';
    end
end

%% Learining parameters
params.epshbias = 1e-1;
params.epsvbias = 1e-1;
params.epsW = 1e-2;
params.phbias = 0.5;
params.pvbias = 0.5;
params.pW = 0.5;
params.decayw = .01;
params.szBatch = 500;
params.sparseness = -1;
params.whitenData = 1;

%% Running parameters
params.iter = 30;
params.verbose = 2;
params.mfIter = 5;
params.saveInterv = 5;
params.useCuda = 0;
params.saveName = 'model.mat';

end
