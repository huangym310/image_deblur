% main_levin.m
%
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone. All rights reserved.

clear all
addpath('lib/');

params.gammaCorrection = false;
params.gamma = 0;
params.iters = 1000;
params.visualize = 1;

dataFolder = 'data/levin/';
resultFolder = 'results/levin/';
if(~exist(resultFolder,'dir'))
    mkdir(resultFolder);
end



lambda = 6e-4;
numIm = 4;
numKer = 8;

for im = 1:numIm
    for ker = 1:numKer
        disp(['im ' num2str(im) ' ker ' num2str(ker)])
        
        name = [dataFolder '/im0' num2str(4 + im) '_flit0' num2str(ker) '.mat'];
        if ~exist(name,'file')
            error(['Levin''s dataset not correctly downloaded. Please ' ...
                   'download it from http://' ...
                   'www.wisdom.weizmann.ac.il/~levina/papers/' ...
                   'LevinEtalCVPR09Data.rar and place all the *.mat ' ...
                   'files in the folder data/levin/.']);
        end
        load(name);
        
        k_ground = f;
        sharp_ground = x;
        blurred = y;
        
        [MK, NK] = size(rot90(k_ground,2));
        tic
        [fe, k] = deblur(blurred, MK, NK, lambda, params);
        t = toc;
        k = rot90(k,2);
        
        % WARNING: The results in the paper are produced using a final
        % non-blind stage with the algorithm presented in the paper "Image and
        % Depth from a Conventional Camera with a Coded Aperture." from
        % Levin et al. You can download the code from:
        % http://groups.csail.mit.edu/graphics/CodedAperture/DeconvolutionCode.html 
        %
        % If you have the code in the MATLAB path, you can use it by
        % uncomment the following line:
        %
        % fe = deconvSps(blurred,k,0.0068);
        
        imshow(fe);
        pause;
        %ssde = comp_upto_shift(fe,sharp_ground);
        
        %resName = [resultFolder 'im0' num2str(4 + im) '_flit0' num2str(ker) '.mat'];
        %save(name,'k', 'fe', 't','ssde');
        
    end
end