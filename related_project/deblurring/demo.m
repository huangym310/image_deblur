function [] = demo()

close all;
addpath bdconv;
I = imread('imgs/lena.png');
I = im2double(I);
I = imresize(I,[150 150]);


img = ones(51,67);
%% 1
k0 = fspecial('gaussian',[9 9],2);
k0 = k0./max(k0(:));
img(1:9,1:9) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);

k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(1:9,11:19) = k;
%% 2
k0 = fspecial('gaussian',[9 9],3);
k0 = k0./max(k0(:));
img(1:9,25:33) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);

k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(1:9,35:43) = k;

%% 3
k0 = zeros(9);
for i=2:8
    k0(10-i,i) = 1;
end
img(1:9,49:57) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);

k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(1:9,59:67) = k;

%% 4
k0 = zeros(9);
for i=2:8
    k0(i,i) = 1;
end
img(15:23,1:9) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);

k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(15:23,11:19) = k;

%% 5
k0 = zeros(9);
k0(2,2:8) = 1;
k0(8,2:8) = 1;
for i=2:8
    k0(10-i,i) = 1;
end
img(15:23,25:33) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(15:23,35:43) = k;
%% 6
k0 = rand(9);
img(15:23,49:57) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(15:23,59:67) = k;

%% 7
k0 = zeros(9);
k0(3:4,6) = 1; 
k0(6:7,4) = 1;
k0(4:6,5) = 1;
img(29:37,1:9) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(29:37,11:19) = k;

%% 8
k0 = zeros(9);
k0(3:7,5) = 1;
img(29:37,25:33) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(29:37,35:43) = k;

%% 9
k0 = zeros(9);
k0(5,3:7) = 1;
img(29:37,49:57) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(29:37,59:67) = k;

%% 10
k0 = zeros(9);
k0(6:8,3) = 1;
k0(4:6,5) = 1;
k0(2:4,7) = 1;
k0(4,5:7) = 1;
k0(6,3:5) = 1;
img(43:51,1:9) = k0;
k0 = k0./sum(sum(k0));
B = conv2(I,k0);
k = estimate_kernel(B,[9 9]);
k = k./max(max(k));
img(43:51,11:19) = k;

figure;
imshow(img);
rmpath bdconv;
end