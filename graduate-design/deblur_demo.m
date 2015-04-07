%clear,clc,close all hidden;
addpath('Levin_data');
%% This is the main script for multi-image deblurring

opt.img_name_set = 'test';
opt.img_format = '.bmp';
opt.img_size = [241 241];
opt.smp_size = [21 21];
opt.ker_size = [15 15];
opt.img_num = 2;

opt.beta = 1.05;
opt.mu_init = 0.00001;
opt.mu_max = 10000;
opt.lambda = 0.1;
opt.max_iter = 200;
opt.epsilon = 1e-10;
opt.sparse = 1;
opt.smooth = 120;
opt.p = 0;

B = prepare_data(opt);
%[B,L] = get_levin_data(opt);

% K = estimate_kernel(B,opt);
% K = estimate_kernel_v2(B,opt);
K = estimate_kernel_v7(B,opt);
S = nonblind_deconv(B,K,opt.lambda);

imshow(S,[]);