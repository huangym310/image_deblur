% In this script we show an example of usage of the core part of the
% deblurring algorithm, as it is presented in the pseudocode of the paper:
%
% D. Perrone and P. Favaro: "Total Variation Blind Deconvolution:
% The Devil is in the Details", IEEE Conference on Computer Vision
% and Pattern Recognition (CVPR), 2014. 
%
% This script does not use the pyramide scheme but only a slow decreasing
% of the parameter lambda after each iteration.
% 
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone. All rights reserved.

clear all
addpath('lib/');

ug = im2double(imread('peppers.png'));
ug = imresize(ug,[111 111]);

kg = fspecial('motion',7);
kg = padarray(kg,[3 0]);
MK = 7;
NK = 7;

f = convn(ug,kg,'valid');

params.iters = 1;
initialLambda = 5e-1;
finalLambda = 1e-4;

u = padarray(f,[floor(MK/2) floor(NK/2)],'replicate');
k = ones(MK,NK)/MK/NK;

lambda = initialLambda;
for i=1:10000
    params.u = u;
    params.k = k;
    [u,k] = blind(f, MK, NK, lambda, params);
    lambda = max(finalLambda,lambda * 0.999);
    
    if mod(i,50) == 0
        figure(1)
        image(uint8([u ug].*255));
        drawnow
        
        figure(2)
        imagesc([k kg])
        drawnow
    end
end
