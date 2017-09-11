function [model, output] = crbm_filtSelection_validation(data, params, utt_num)
% Pass the validation data spectrograms through convolutional restricted Boltzmann machine 
%   with the weight being fixed as learned filters.
%
%   [model output] = crbm_filtSelection_validation(data, params, utt_num)
%
%   data should be a structure, containing:
%       data.x      The input sampled from spectrogram. This matrix is 4-D.  The first two dimensions
%                   define the size of the input sampled from spectrogram (third dimension is color channel, which is one in our case),
%                   and the last dimension indexes through the batch of
%                   images. I.e. the four dimensions are: height, width,
%                   channels (1 for grayscale, 3 for RGB), and number of
%                   sampled inputs.
%
%   Written by: Purvi Agrawal
%   January, 2017
%   LEAP Lab, IISc
%
%   Original code can be found at : http://qipeng.me/software/convolutional-rbm.html
%   https://github.com/qipeng/convolutionalRBM.m. For research / personal
%   purposes only.

if params.verbose > 0,
    fprintf('Starting training CRBM with the following parameters:\n');
    disp(params);
    fprintf('Initializing parameters...');
end

useCuda = params.useCuda;

if isfield(params, 'method'),
    if strcmp(params.method, 'CD'),
        method = 1; % Contrastive Divergence
    elseif strcmp(params.method, 'PCD'),
        method = 2; % Persistent Contrastive Divergence
    end
else
    method = 1;     % use Contrastive Divergence as default
end

%% initialization
N = size(data, 4);
Nfilters = params.nmap;
Wfilter = params.szFilter;
Hfilter=Wfilter(1); Wfilter=Wfilter(2);

p = params.szPool;
H = size(data, 1);
W = size(data, 2);
colors = size(data, 3);
Hhidden = H - Hfilter + 1;
Whidden = W - Wfilter + 1;
Hpool = floor(Hhidden / p);
Wpool = floor(Whidden / p);
param_iter = params.iter;
% param_szBatch = params.szBatch;
param_szBatch = size(data,4); %*******************
output_enabled = nargout > 1;

%% Load the filter (rate or scale) for which the hidden activations are to be evaluated

dataFile = 'data/model_rate/model_clean_mel_rate_15tap_1.mat';
% dataFile = 'data/model_rate/model_clean_mel_rate_15tap_2.mat';
% dataFile = 'data/model_rate/model_clean_mel_rate_15tap_3.mat';
% dataFile = 'data/model_scale/model_clean_mel_scale_8tap_1.mat';
% dataFile = 'data/model_scale/model_clean_mel_scale_8tap_2.mat';
% dataFile = 'data/model_scale/model_clean_mel_scale_8tap_3.mat';

w1 = load(dataFile);
model.W = w1.model.W;
model.vbias = w1.model.vbias;
model.hbias = w1.model.hbias;
model.sigma = w1.model.sigma;

if output_enabled,
    output.x = zeros(Hpool, Wpool, Nfilters, N);
end

total_batches = ceil(N / param_szBatch);

if params.verbose > 0,
    fprintf('Completed.\n');
end

if ~isfield(model,'iter')
    model.iter = 0;
end

%% Computing one mean and std of full data
% mu=mean(data.x(:)); std_dev=std(data.x(:));
% % save('data/model_650_50_lrchange/data_2_mean_std','mu', 'std_dev');
% 
data = bsxfun(@rdivide, bsxfun(@minus, data, mean(data(:))), std(data(:)));

%% COmputing mean and std col. wise
% mu=mean(data); std_dev=std(data);
% mu_rev=repmat(mu,[size(data,1),1,1,1]); 
% std_dev_rev=repmat(std_dev,[size(data,1),1,1,1]);
% data=(data-mu_rev)./std_dev_rev;

%%
if method == 2,
    phantom = randn(H, W, colors, N);
end
error_iter=zeros(1,param_iter);

for iter = model.iter+1:param_iter,
    % shuffle data
    batch_idx = randperm(N);
    
    if params.verbose > 0,
        fprintf('Iteration %d\n', iter);
        if params.verbose > 1,
            fprintf('Batch progress (%d total): ', total_batches);
        end
    end
    
    errsum = 0;
    
    if (iter > 5),
        params.pW = .9;
        params.pvbias = 0;
        params.phbias = 0;
    end
    
    for batch = 1:total_batches,
        batchdata = data(:,:,:,batch_idx((batch - 1) * param_szBatch + 1 : ...
            min(batch * param_szBatch, N)));
        if method == 2,
            phantomdata = phantom(:,:,:,((batch - 1) * param_szBatch + 1 : ...
                min(batch * param_szBatch, N)));
        end
        recon = batchdata;
        
        %% positive phase

        %% hidden update
        
        model_W = model.W;
        model_hbias = model.hbias;
        model_vbias = model.vbias;
        
        poshidacts = convs(recon, model_W, useCuda);

         [poshidprobs, pospoolprobs, poshidstates] = poolHidden(poshidacts / model.sigma, model_hbias / model.sigma, p, useCuda);
%        [poshidprobs, pospoolprobs, poshidstates] = poolH_matlab(poshidacts / model.sigma, model_hbias / model.sigma, p);
     
        %% negative phase
        
        %% reconstruct data from hidden variables

        if method == 1,
            recon = conve(poshidprobs, model_W, useCuda);
        elseif method == 2,
            recon = phantomdata;
        end
        
        recon = bsxfun(@plus, recon, reshape(model_vbias, [1 1 colors]));

        if (params.sparseness > 0),
            recon = recon + model.sigma * randn(size(recon));
        end
        
        %% mean field hidden update
        
        neghidacts = convs(recon, model_W, useCuda);
        neghidprobs = poolHidden(neghidacts / model.sigma, model_hbias / model.sigma, p, useCuda);
            
        if (params.verbose > 1),
            fprintf('.');
            err = batchdata - recon;
            errsum = errsum + sum(err(:).^2);             
        end
    end
    
    if (params.verbose > 1),
        fprintf('\n\terror:%f', errsum);
        error_iter(iter)=errsum;
    end
    error_iter(iter)=errsum;
    
    if (iter == param_iter),

            model.iter = iter;
            dest= strcat('data/model_hid_prob','/','clean_mel_15tap_rate1'); % make a destination folder corresponding to the loaded filter to save hidden activation probabilities for each validation data spectrogram
            mkdir(dest)
            save(strcat(dest,'/','model_clean_mel_rate1_utt_',num2str(utt_num)), 'model', 'iter','error_iter', 'errsum','poshidprobs','neghidprobs');
    end   
end
