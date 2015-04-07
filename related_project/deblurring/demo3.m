function [] = demo3()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B = imread('imgs/15.jpg');
B = rgb2gray(B);
B = im2double(B);
B = B(700:1050,730:1080);
imwrite(B,'results/R1.png');

ksize = [31 31];rate = 0.75; lambda = 0.0015;alpha = 27000;
addpath bdconv;
[k,Id] = bldconv_sp(B,ksize,rate,alpha,lambda);
rmpath bdconv;
    
imwrite(Id,'results/R1_deblur.png');
imwrite(k./max(k(:)),'results/R1_kernel.png');

end

