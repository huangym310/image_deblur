% main_noisy.m
%
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone. All rights reserved.

clear all
addpath('lib/');
name = 'noisy';

%path variables
path = 'data/';
blurredPath = 'noisy.jpg';
resultsPath = 'results/';

if(~exist(resultsPath,'dir'))
    mkdir(resultsPath)
end

%% -----------------------------------------------------------------------%
%                                  Input                                  %
% ------------------------------------------------------------------------%
g = imread([path blurredPath]);

%% -----------------------------------------------------------------------%
%                                Parameters                               %
% ------------------------------------------------------------------------%
MK = 31;
NK = 45;
lambda = 1e-3;

%% -----------------------------------------------------------------------%
%                                  Deblur                                 %
% ------------------------------------------------------------------------%
[u, k] = deblur(g,MK,NK,lambda);

u(u<0) = 0;
u(u>1) = 1;
imwrite(u,[resultsPath name 'out.png']);
imwrite(imresize(k./max(k(:)),5,'nearest'),[resultsPath name '_kernel.png']);