function [u, k] = coarseToFine(f, MK, NK, blind_params, params)
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone.  All rights reserved.

u = padarray(f,[floor(MK/2) floor(NK/2)],'replicate');
k = ones(MK,NK)/MK/NK;

[fp, Mp, Np, MKp, NKp, lambdas, scales] = buildPyramid(f, MK, NK,...
                params.finalLambda, params.lambdaMultiplier, ...
                params.interpolationMethod, params.kernelSizeMultiplier,...
                params.maxLambda);

            

% Multiscale Processing
for scale=scales:-1:1
    Ms = Mp{scale,1};
    Ns = Np{scale,1};
    
    MKs = MKp{scale,1};
    NKs = NKp{scale,1}; 

    u = imresize(u, [(Ms + MKs - 1) (Ns + NKs - 1)], 'Method',...
        params.interpolationMethod);

    k = imresize(k, [MKs NKs], 'Method', params.interpolationMethod);
    k = k.*(k > 0);
    k = k./sum(k(:));
    
    fs = fp{scale,1};
    
    lambda = lambdas(scale);
    
    blind_params.u = u;
    blind_params.k = k;
    
    if (blind_params.visualize) 
        disp(['scale: ' num2str(scale) ' lambda: ' num2str(lambda) ' MKs: '...
                num2str(MKs) ' NKs: ' num2str(NKs) ' iters: ' num2str(blind_params.iters)])
    end
    
    [u, k] = blind(fs, MKs, NKs, lambda, blind_params);
    u = dec(fs, k, lambda, blind_params);
end
