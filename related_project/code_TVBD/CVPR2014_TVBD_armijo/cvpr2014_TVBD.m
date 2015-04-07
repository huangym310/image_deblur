function [L,K] = cvpr2014_TVBD(B,ks,lambda)
if ~exist('lambda','var')
    lambda = 6e-4;
end
if size(B,3) > 1
    B_gray = rgb2gray(B);
else 
    B_gray = B;
end
if max(B_gray(:))>2
    B_gray = double(B_gray)/255;
    B = double(B)/255;
else 
    B_gray = double(B_gray);
    B = double(B);
end
n = size(B_gray,3);
%% paramsters setting
% adaptive lambda with noise
sigma = 0;% for noise level
filt_size = 5;
filt = ones(filt_size)/filt_size^2;
for ii = 1:n
    TVBi = TV(B(:,:,ii));
    Bi_smooth = conv2(B(:,:,ii),filt,'same');
    TVBi_smooth = TV(Bi_smooth);   
    sigma = sigma+TVBi/TVBi_smooth/n;    
end
% lambda in [3e-4,6e-4],with noise [1e-3,3e-3]
lambda = lambda*sigma;

params = get_params_multi(lambda);
%% multi-scale framework to estimate blur kernel
K = estimate_kernel_multi_scale(B_gray,ks,params);
%% final deblur
L = non_blind_decov(B,K,params.lambda,500);
hkr = floor(ks(1)/2);% kernel size
hkc = floor(ks(2)/2);% kernel size
L = L(hkr+1:end-hkr,hkc+1:end-hkc,:);
L = min(max(L,0),1);%[0,1]


%end

