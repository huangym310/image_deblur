function [] = demo1()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ksize = [13 13];rate = 0.75; lambda = 0.001;alpha = 1050;
B = imread('imgs/I1.png');
B = im2double(B);
   
addpath bdconv;
[k,Id] = bldconv_sp(B,ksize,rate,alpha,lambda);
rmpath bdconv;
    
imwrite(Id,'results/I1_deblur.png');
imwrite(k./max(k(:)),'results/I1_kernel.png');

end

