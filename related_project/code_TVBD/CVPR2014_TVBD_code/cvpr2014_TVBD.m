function [L,K] = cvpr2014_TVBD(B,ks,lambda)
if ~exist('lambda','var')
    lambda = 6e-4;
end

if size(B,3) > 1
    B_gray = rgb2gray(B);
else 
    B_gray = B;
end
B_gray = im2double(B_gray);
%% paramsters setting
% adaptive lambda with noise
Bx_gray = B_gray(:,[2:end,end])-B_gray;
By_gray = B_gray([2:end,end],:)-B_gray;
filt = ones(5)/5^2;% mean filter denoising
B_gray_smooth = conv2(B_gray,filt,'same');
Bx_gray_smooth = B_gray_smooth(:,[2:end,end])-B_gray_smooth;
By_gray_smooth = B_gray_smooth([2:end,end],:)-B_gray_smooth;
Bg_gray = sqrt(Bx_gray.^2+By_gray.^2);
Bg_gray_smooth = sqrt(Bx_gray_smooth.^2+By_gray_smooth.^2);
% lambda in [3e-4,6e-4],with noise [1e-3,3e-3]
lambda = lambda*sum(Bg_gray(:))/sum(Bg_gray_smooth(:));

params = get_params_multi(lambda);
%% multi-scale framework to estimate blur kernel
K = estimate_kernel_multi_scale(B_gray,ks,params);
%% final deblur
L = non_blind_decov(B,K,params.lambda,params.nloops);
kr = ks(1);% kernel size
kc = ks(2);
L = L(floor(kr/2)+1:end-floor(kr/2),floor(kc/2)+1:end-floor(kc/2),:);
L = min(max(L,0),1);%[0,1]


end

