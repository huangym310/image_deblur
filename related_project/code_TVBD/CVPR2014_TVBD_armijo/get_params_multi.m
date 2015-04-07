function params = get_params_multi(lambda)
%% initialization for multi-scales
params.ks_min = 3;% the coarsest kernel size
params.scale_factor = 1.1;% the scale factor of pyramid
%% for Alternative Minization
params.nloops = 1000;% the loops for each level
params.interpMethod = 'bicubic';%'nearest', 'bilinear', 'bicubic'

%% for update image
params.lambda = lambda;% weight for L2 norm of reweighted image prior
params.lambda_max = 1.1e-1;% weight for L2 norm of blur kernel
params.lambda_step = 1;
%% for update kernel
params.kernel_threshold = 0/100;% of the max kernel value
params.initializeMethod = 'ver';% 'ver','hor','uniform','delta'
params.gamma = 0;
end

