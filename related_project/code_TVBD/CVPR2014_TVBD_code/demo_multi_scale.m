clear,clc,close all hidden;
load('im05_flit01.mat');
% K_gt = f;% ground true kernel
% L_gt = x;% ground true image
% %% input
% B = y;%conv2(L_gt,K_gt,'same');% blurry image
K_gt = zeros(19,19);
K_gt(4:13,10) = 0.05;
K_gt(4,11:15) = 0.05;
K_gt(13,5:9) = 0.05;
L_gt = x;
B = conv2(L_gt,K_gt,'valid');
ks = [19,19];% the blur kernel size

[L,K] = cvpr2014_TVBD(B,ks);

figure,subplot(1,2,1);imshow(mat2gray(K_gt));title('Ground Truth');
subplot(1,2,2);imshow(mat2gray(K));title('Estimated Kernel');

figure,subplot(1,3,1);imshow(mat2gray(L_gt));title('Sharp Image');
subplot(1,3,2);imshow(mat2gray(B));title('Blurry Image'); 
subplot(1,3,3);imshow(mat2gray(L));title('Deblurred Image');               