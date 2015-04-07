function [u, k] = deblur(f,MK,NK,lambda,params)
% [u, k] = deblur(f,MK,NK,lambda) implements the blind
% deconvolution algorithm presented in the paper:
%
% D. Perrone and P. Favaro: "Total Variation Blind Deconvolution:
% The Devil is in the Details", IEEE Conference on Computer Vision
% and Pattern Recognition (CVPR), 2014. 
%
% Where f is the input blurry image, MK and NK are the support of
% the unknown PSF, and lambda weights the total variation
% regularization. In output the function gives the sharp image u
% and the PSF k. 
% A typical choice for lambda is between 3e-4 and 6e-4. For severly
% noisy images it might be necessary to use values between 1e-3 and
% 3e-3.  
%
% [u, k] = deblur(f,MK,NK,lambda, params) allows the user to
% specify in the struct params additional input parameters, as:
%
% params.gammaCorrection:: [bool]
%        If true, performs gamma correction before the execution of
%        the deblurring algorithm.
%
% params.gamma:: [double]
%        Value for gamma correction.
%
% params.iters:: [uint(1000)]
%        Number of innermost iterations for the gradient descent
%        algorithm.
%
% params.visualize:: [0-1]
%        If set to 1 visualize diagnostics and current estimates.
%
% Example: 
% ug = im2double(imread('peppers.png'));
% kg = fspecial('motion',15);
% f = convn(ug,kg,'valid');
% [u, k] = deblur(f,11,15,3e-4);
% 
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone. All rights reserved.


if (~exist('params','var'))
    params.gammaCorrection = false;
    params.gamma = 0;
    params.iters = 1000;
    params.visualize = 1;
else
    if ~isfield(params,'iters')
        params.iters = 1000;
    end
    if ~isfield(params,'gammaCorrection')
        params.gammaCorrection = false;
    end
    if ~isfield(params,'gamma')
        params.gamma = 0;
    end
    if ~isfield(params,'visualize')
        params.visualize = 0;
    end
end

%% -------------------------------------------------------------------%
%                           Initialization                            %
% --------------------------------------------------------------------%
% Convert B to double [0,1]
f = im2double(f);
    
[M, N, C] = size(f);

if(C == 1)
    % make sure that the image size is an odd number (on both axes)
    if (mod(N,2) == 0)
        f = f(1:end-1,:);
    end
    if (mod(M,2) == 0)
        f = f(:,1:end-1);
    end
else
    % make sure that the image size is an odd number (on both axes)
    if (mod(N,2) == 0)
        f = f(1:end-1,:,:);
    end
    if (mod(M,2) == 0)
        f = f(:,1:end-1,:);
    end
end


%% -------------------------------------------------------------------%
%                            Preprocessing                            %
% --------------------------------------------------------------------%
if (params.gammaCorrection)
    f = gammaCorrection(f,gamma);
end

%% -------------------------------------------------------------------%
%                            PSF Estimation                           %
% --------------------------------------------------------------------%
% coarse to fine parameters
ctf_params.lambdaMultiplier = 1.9;
ctf_params.maxLambda = 1.1e-1;
ctf_params.finalLambda = lambda;
ctf_params.kernelSizeMultiplier = 1.1;
ctf_params.interpolationMethod = 'bicubic';

[u, k] = coarseToFine(f, MK, NK, params, ctf_params);
