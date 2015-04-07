function [] = demo2()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ksize = [13 13];rate = 0.75; lambda = 0.002;alpha = 5500;
B = imread('imgs/I2.png');
B = im2double(B);
   
addpath bdconv;
[k,Id] = bldconv_sp(B,ksize,rate,alpha,lambda);
rmpath bdconv;
    
imwrite(Id,'results/I2_deblur.png');
imwrite(k./max(k(:)),'results/I2_kernel.png');

end

