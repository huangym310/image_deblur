% main_postcard.m
%
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone. All rights reserved.

clear sll
addpath('lib/');
name = 'postcard';

%path variables
path = 'data/';
blurredPath = 'postcard.png';
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
MK = 75;
NK = 35;
lambda = 3.5e-4;

%% -----------------------------------------------------------------------%
%                                  Deblur                                 %
% ------------------------------------------------------------------------%
[u, k] = deblur(g,MK,NK,lambda);

u(u<0) = 0;
u(u>1) = 1;
imwrite(u,[resultsPath name 'out.png']);
imwrite(imresize(k./max(k(:)),5,'nearest'),[resultsPath name '_kernel.png']);

