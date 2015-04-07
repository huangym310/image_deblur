clear,clc,close all hidden;
load('im05_flit01.mat');
K_gt = f;% ground true kernel
L_gt = x;% ground true image
%% input
B = y;%conv2(L_gt,K_gt,'same');% blurry image
ks = [19,19];% the blur kernel size

[L,K] = cvpr2014_TVBD(B,ks);
% cvpr2014_TVBD;

figure,subplot(1,2,1);imshow(mat2gray(K_gt));title('Ground Truth');
subplot(1,2,2);imshow(mat2gray(K));title('Estimated Kernel');

figure,subplot(1,3,1);imshow(mat2gray(L_gt));title('Sharp Image');
subplot(1,3,2);imshow(mat2gray(B));title('Blurry Image'); 
subplot(1,3,3);imshow(mat2gray(L));title('Deblurred Image');   