function [model, output] = trainCRBM_Rate(data, params)
% trainCRBM_Rate  Trains a convolutional restricted Boltzmann machine 
%   with the specified parameters.
%
%   [model output] = trainCRBM_Rate(data, params)
%
%   data should be a structure, containing:
%   data.x      The input sampled from spectrogram. This matrix is 4-D.  The first two dimensions
%               define the size of the input sampled from spectrogram for respective filter learning (third dimension is color channel, which is one in our case),
%               and the last dimension indexes through the batch of
%               images. I.e. the four dimensions are: height, width,
%               color channels (1 for grayscale, 3 for RGB), and number of
%               sampled inputs.
%
%
%   Written by: Purvi Agrawal, 
%   September, 2016
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
N = size(data.x, 4);
Nfilters = params.nmap;
Wfilter = params.szFilter;
Hfilter=Wfilter(1); Wfilter=Wfilter(2);

p = params.szPool;
H = size(data.x, 1);
W = size(data.x, 2);
colors = size(data.x, 3);
Hhidden = H - Hfilter + 1;
Whidden = W - Wfilter + 1;
Hpool = floor(Hhidden / p);
Wpool = floor(Whidden / p);
param_iter = params.iter;
param_szBatch = params.szBatch;
output_enabled = nargout > 1;

%vmasNfilters = conve(ones(nh), ones(m), useCuda);

hinit = 0;

if params.sparseness > 0,
    hinit = -.2;
end


if exist('model.mat') && ~isempty('model.mat'),
    load('model.mat')
    model = model;
    if (~isfield(model,'W')), 
        model.W = 0.01 * randn(Hfilter, Wfilter, colors, Nfilters);
    else
        if (size(model.W) ~= [Hfilter Wfilter colors Nfilters]), error('Incompatible input model.'); end
    end
    if (~isfield(model,'vbias')), model.vbias = zeros(1, colors);end
    if (~isfield(model,'hbias')), model.hbias = ones(1, Nfilters) * hinit;end
    if (~isfield(model,'sigma')),
        if (params.sparseness > 0)
            model.sigma = 0.1;
        else
            model.sigma = 1;    
        end
    end
else
    %initializing filter
    model.W = 0.01 * randn(Hfilter, Wfilter, colors, Nfilters);
    
    model.vbias = zeros(1, colors);
    model.hbias = ones(1, Nfilters) * hinit;
    if (params.sparseness > 0)
        model.sigma = 0.1;
    else
        model.sigma = 1;    
    end
end

dW = 0;
dvbias = 0;
dhbias = 0;

pW = params.pW;
pvbias = params.pvbias;
phbias = params.phbias;

if output_enabled,
    output.x = zeros(Hpool, Wpool, Nfilters, N);
end

total_batches = ceil(N / param_szBatch);

if params.verbose > 0,
    fprintf('Completed.\n');
end

hidq = params.sparseness;
lambdaq = 0.9;

if ~isfield(model,'iter')
    model.iter = 0;
end

%% Computing one mean and std of full data
% mu=mean(data.x(:)); std_dev=std(data.x(:));
% save('data/model_650_50_lrchange/data_2_mean_std','mu', 'std_dev');

% data.x = bsxfun(@rdivide, bsxfun(@minus, data.x, mean(data.x(:))), std(data.x(:)));
% data.x = bsxfun(@minus, data.x, mean(data.x(:)));
%% Computing mean and std col. wise
mu=mean(data.x);
% std_dev=std(data.x);
mu_rev=repmat(mu,[size(data.x,1),1,1,1]); 
% std_dev_rev=repmat(std_dev,[size(data.x,1),1,1,1]);
% data.x=(data.x-mu_rev)./std_dev_rev;
data.x=(data.x-mu_rev);

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
    
    hidact = zeros(1, Nfilters);
    errsum = 0;
    
    if (iter > 5),
        params.pW = .9;
        params.pvbias = 0;
        params.phbias = 0;
    end
    
    for batch = 1:total_batches,
        batchdata = data.x(:,:,:,batch_idx((batch - 1) * param_szBatch + 1 : ...
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
            if (params.verbose > 4),
                %% visualize data, reconstruction, and filters (still experimental)
                figure(10);
                for i = 1:1,subplot(3,1,i+1);imagesc(model.W(:,:,:,i));axis image off;end;colormap default;drawnow;
                subplot(2,2,1);imagesc(batchdata(:,:,1));colormap default;axis off;title('data (ZCA''d)');
                subplot(2,2,2);imagesc(recon(:,:,1));colormap default;axis off;title('reconstruction');
                drawnow;
            end
        end
        
        %% contrast divergence update on params
        
        if (params.sparseness > 0),  % for sparsity, we need to consider probabilities after pooling (hence, pospoolprobs): 
            hidact = hidact + reshape(sum(sum(sum(pospoolprobs, 4), 2), 1), [1 Nfilters]);  % with sparsity constraint, only hidden bias is changed in our case
        % (here, accumulating hidden activation probabilities of all the previous batches till current batch)
        else
            dhbias = phbias * dhbias + ...
                reshape((sum(sum(sum(poshidprobs, 4), 2), 1) - sum(sum(sum(neghidprobs, 4), 2), 1))...
                / Whidden / Hhidden / param_szBatch, [1 Nfilters]);
        end
        
        dvbias = pvbias * dvbias + ...
            reshape((sum(sum(sum(batchdata, 4), 2), 1) - sum(sum(sum(recon, 4), 2), 1))...
            / H / W / param_szBatch, [1 colors]);
       ddw = convs4(batchdata(Hfilter:H-Hfilter+1,Wfilter:W-Wfilter+1,:,:), poshidprobs(Hfilter:Hhidden-Hfilter+1,Wfilter:Whidden-Wfilter+1,:,:), useCuda) ...
            - convs4(    recon(Hfilter:H-Hfilter+1,Wfilter:W-Wfilter+1,:,:), neghidprobs(Hfilter:Hhidden-Hfilter+1,Wfilter:Whidden-Wfilter+1,:,:), useCuda);
      % here, using only one step of Gibbs sampling for training RBM
        dW = pW * dW + ddw / (Hhidden - 2 * Hfilter + 2) / (Whidden - 2 * Wfilter + 2) / param_szBatch;
        
        %%*********saving previous model parameters ***
        model.vbias_prev = model.vbias ;
        model.hbias_prev = model.hbias ;
        model.W_prev = model.W ;
        % ***************************
        
        model.vbias = model.vbias + params.epsvbias * dvbias;
        if params.sparseness <= 0, % i.e. no sparsity constraint:
            model.hbias = model.hbias + params.epshbias * dhbias; 
        end
        model.W = model.W + params.epsW * (dW  - params.decayw * model.W);
        
    end
    
    if (params.verbose > 1),
        fprintf('\n\terror:%f', errsum);
        error_iter(iter)=errsum;
    end
    error_iter(iter)=errsum;
    %%%%%%%% added to compensate for increased error at current iteration
    if iter > 1
        if error_iter(iter) > error_iter(iter - 1 )
             params.epsW = params.epsW/2 ;  % reducing the learning rate according to the error
             params.epsvbias = params.epsvbias/2 ;
             params.epshbias = params.epshbias/2 ;
             model.W = model.W_prev ;
             model.hbias = model.hbias_prev;
             model.vbias = model.vbias_prev ;
            disp(['******** Error Increased ' num2str(params.epsW) ' **********' ]);
        end
    end
    %%%%%%%%%
    if params.sparseness > 0,
        hidact = hidact / Hhidden / Whidden / N;  %% 'Mean' activation probability of the hidden unit till the current mini-batch
        
        hidq = hidq * lambdaq + hidact * (1 - lambdaq); % lambdaq = 0.9 is the decay rate, it remains fixed. hidq is the new sparseness parameter getting updated according to
        % the mean hidden activations (hidact), lambdaq and the old sparseness
        % parameter (hidq). New hidq is estimated by using an exponentially
        % decaying avg. of the mean probability that a unit is active in
        % each mini-batch.
        
        dhbias = phbias * dhbias + ((params.sparseness) - (hidq)); % change in the sparsity parameter is considered for updating hidden bias.
        
        model.hbias = model.hbias + params.epshbias * dhbias;
        if params.verbose > 0,
            if (params.verbose > 1),
                fprintf('\tsigma:%f', model.sigma);
            end
            fprintf('\n\tsparseness: %f\thidbias: %f\n', sum(hidact) / Nfilters, sum(model.hbias) / Nfilters);
        end
        if (model.sigma > 0.01),
            model.sigma = model.sigma * 0.99;
        end
    end
    
    if ~rem(iter, param_iter),

            model.iter = iter;
            % For 1st Rate filter learning (otherwise comment)
            save('data/model_rate/model_clean_mel_rate_15tap_1', 'model', 'iter','error_iter');
            
            % For 2nd Rate filter learning (otherwise comment)
            % save('data/model_rate/model_clean_mel_rate_15tap_2', 'model', 'iter','error_iter');
            
            % For 3rd Rate filter learning (otherwise comment)
            % save('data/model_rate/model_clean_mel_rate_15tap_3', 'model', 'iter','error_iter');
            
            if params.verbose > 1,
                fprintf('Model saved at iteration %d\n', iter);
            end
%         end
    end
    
   
end

% figure;
% plot(error_iter); title('Error at each iteration');
