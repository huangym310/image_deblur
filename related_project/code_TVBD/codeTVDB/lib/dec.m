function u = dec(f, k, lambda, params)
% u = dec(f, MK, NK, lambda) implements the non-blind
% deconvolution gradient decent algorithm presented in the paper:
%
% D. Perrone and P. Favaro: "Total Variation Blind Deconvolution:
% The Devil is in the Details", IEEE Conference on Computer Vision
% and Pattern Recognition (CVPR), 2014. 
%
% Where f is the input blurry image,k is the PSF, and lambda
% weights the total variation regularization. In output the
% function gives the sharp image u.  
% A typical choice for lambda is between 3e-4 and 6e-4. For severly
% noisy images it might be necessary to use values between 1e-3 and
% 3e-3.  
%
% u = deblur(f, k, lambda, params) allows the user to
% specify in the struct params additional input parameters, as:
%
% params.u:: [(f)]
%        Initial sharp image.
%
% params.k:: [(ones(MK,NK)/MK/NK)]
%        Initial PSF.
%
% params.iters:: [uint(1000)]
%        Number of innermost iterations for the gradient descent
%        algorithm.
%
% params.visualize:: [0-1]
%        If set to 1 visualize diagnostics and current estimates.
% 
% Authors: Daniele Perrone, perrone@iam.unibe.ch
%          Paolo Favaro, paolo.favaro@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone and Paolo Favaro.  All rights
% reserved. 

% get image dimension
[M, N, C] = size(f);
[MK, NK] = size(k);

if ~exist('params','var')
    params.iters = 1000;
    params.u = padarray(f,[floor(MK/2) floor(NK/2)],'replicate');
    params.visualize = 0;
else
    if ~isfield(params,'iters')
        params.iters = 1000;
    end
    if ~isfield(params,'u')
        params.u = padarray(f,[floor(MK/2) floor(NK/2)],'replicate');
    end
    if ~isfield(params,'visualize')
        params.visualize = 0;
    end
end

if ~exist('lambda','var')
    lambda = 3e-4;
end

% size of sharp image
MU = M + MK - 1;
NU = N + NK - 1;
gradudata = zeros(MU,NU,C);
u = params.u;

for  it = 1:params.iters
    for c=1:C
        if numel(k)>41^2
            gradudata(:,:,c) = conv2fft(conv2fft(u(:,:,c), k ,'valid') ...
                                        -  f(:,:,c), rot90(k,2), 'full');
        else
            gradudata(:,:,c) = conv2(conv2(u(:,:,c), k, 'valid') - ...
                                     f(:,:,c), rot90(k,2), 'full');
        end
    end
    
    gradu = (gradudata - lambda*gradTVcc(u));
    sf = 5e-3*max(u(:))/max(1e-31,max(max(abs(gradu(:)))));
    u   = u - sf*gradu;
    
    if params.visualize
        if mod(it,50)==0
            figure(1);
            uv = u;
            uv(uv > 1) = 1;
            uv(uv <0 ) = 0;
            fv=  f;
            fv(fv > 1) = 1;
            fv(fv < 0) = 0;
            fv = padarray(fv,[floor(MK/2) floor(NK/2)]);
            imagesc([uv fv]);
            drawnow;
        end
    end
end
