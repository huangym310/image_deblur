function [] = demo4()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B = imread('imgs/19.jpg');
B = rgb2gray(B);
B = im2double(B);
B = B(600:900,1200:1500);
imwrite(B,'results/R2.png');
%% 
ksize = [17 17];rate = 0.75; lambda = 0.0012;alpha = 70;

addpath bdconv;
[k,Id] = bldconv_sp(B,ksize,rate,alpha,lambda);
rmpath bdconv;
    
imwrite(Id,'results/R2_deblur.png');
imwrite(k./max(k(:)),'results/R2_kernel.png');

end


