function [] = demo5()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B = imread('imgs/20.jpg');
B = rgb2gray(B);
B = im2double(B);
B = B(50:1050,800:1900);
B = imresize(B,[200 200]);
imwrite(B,'results/R3.png');
%% 
ksize = [17 17];rate = 0.75; lambda = 0.001;alpha = 4000;

addpath bdconv;
[k,Id] = bldconv_sp(B,ksize,rate,alpha,lambda);
rmpath bdconv;
    
imwrite(Id,'results/R3_deblur.png');
imwrite(k./max(k(:)),'results/R3_kernel.png');

end


