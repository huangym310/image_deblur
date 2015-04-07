function S = nonblind_deconv(B,kernel,lambda)
%%
% Image restoration with L0 prior
% The objective function: 
% S^* = argmin \sum_i^n 1/n ||I*k_i - B_i||^2 + lambda |\nabla I|_0
% qual to
% argmin \sum_i^n 1/n 1/lambda ||I*k_i - B_i||^2 + |\nabla I|_0
%% Input:
% @Im: Blurred image
% @kernel: blur kernel
% @lambda: weight for the L0 prior
% @kappa: Update ratio in the ADM
%% Output:
% @S: Latent image
%
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang,
%        Deblurring Text Images via L0-Regularized Intensity and Gradient
%        Prior, CVPR, 2014. 
%   [2] Li Xu, Cewu Lu, Yi Xu, and Jiaya Jia. Image smoothing via l0 gradient minimization.
%        ACM Trans. Graph., 30(6):174, 2011.
%
%   Author: Jinshan Pan (sdluran@gmail.com)
%   Date  : 05/18/2014

%% pad image
H = size(B,1);    W = size(B,2);
B = wrap_boundary(B,opt_fft_size([H+size(kernel,1),W+size(kernel,2)]-1));
%% pre compuate
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(B);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
FB = B;
for ii = 1:D
    FB(:,:,ii) = fft2(B(:,:,ii));
end
%% 
n = D;
S = mean(B,3);
Den_KER = zeros(sizeI2D);
Normin1 = zeros(sizeI2D);
for ii = 1:D
    this_KER = psf2otf(kernel(:,:,ii),sizeI2D);
    Den_KER = Den_KER+abs(this_KER).^2/n;    
    Normin1 = Normin1+conj(this_KER).*FB(:,:,ii)/n;
end
%% 
betamax = 1e5;
beta = min(1,2*lambda);
beta_step = 2;
while beta < betamax
    %h-v
    Denormin  = Den_KER + beta*Denormin2;
    h = [diff(S,1,2), S(:,1,:)-S(:,end,:)];
    v = [diff(S,1,1); S(1,:,:)-S(end,:,:)];
    t = (h.^2+v.^2)<lambda/beta;    
    h(t)=0; 
    v(t)=0;
    % get S
    Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
    Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
    FS = (Normin1 + beta*fft2(Normin2))./Denormin;
    S = real(ifft2(FS));   
    % update beta
    beta = beta*beta_step;
end
S = S(1:H,1:W,:);

end

% function x = nonblind_deconv(I,filt1,we,max_it)
% %note: size(filt1) is expected to be odd in both dimensions 
% if (~exist('max_it','var'))
%    max_it  =200;
% end
% 
% [n,m]=size(I);
% hfs1_x1=floor((size(filt1,2)-1)/2);
% hfs1_x2=ceil((size(filt1,2)-1)/2);
% hfs1_y1=floor((size(filt1,1)-1)/2);
% hfs1_y2=ceil((size(filt1,1)-1)/2);
% 
% hfs_x1=hfs1_x1;
% hfs_x2=hfs1_x2;
% hfs_y1=hfs1_y1;
% hfs_y2=hfs1_y2;
% 
% m=m+hfs_x1+hfs_x2;
% n=n+hfs_y1+hfs_y2;
% 
% 
% tI=I;
% I=zeros(n,m);
% I(hfs_y1+1:n-hfs_y2,hfs_x1+1:m-hfs_x2)=tI; 
% x=I;
% 
% dxf=[1 -1];
% dyf=[1;-1];
% dyyf=[-1; 2; -1];
% dxxf=[-1, 2, -1];
% dxyf=[-1 1;1 -1];
% 
% weight_x=ones(n,m-1);
% weight_y=ones(n-1,m);
% weight_xx=ones(n,m-2);
% weight_yy=ones(n-2,m);
% weight_xy=ones(n-1,m-1);
% 
% x = deconvL2_w(x(hfs_y1+1:n-hfs_y2,hfs_x1+1:m-hfs_x2),filt1,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy);
% 
% w0=0.1;
% exp_a=0.8;
% thr_e=0.01; 
% for t=1:2
%   dy=conv2(x,rot90(dyf,2),'valid');
%   dx=conv2(x,rot90(dxf,2),'valid');
%   dyy=conv2(x,rot90(dyyf,2),'valid');
%   dxx=conv2(x,rot90(dxxf,2),'valid');
%   dxy=conv2(x,rot90(dxyf,2),'valid');
%   
%   weight_x=w0*max(abs(dx),thr_e).^(exp_a-2); 
%   weight_y=w0*max(abs(dy),thr_e).^(exp_a-2);
%   weight_xx=0.25*w0*max(abs(dxx),thr_e).^(exp_a-2); 
%   weight_yy=0.25*w0*max(abs(dyy),thr_e).^(exp_a-2);
%   weight_xy=0.25*w0*max(abs(dxy),thr_e).^(exp_a-2);
%   
%   x = deconvL2_w(I(hfs_y1+1:n-hfs_y2,hfs_x1+1:m-hfs_x2),filt1,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy);
% 
% end
% x=x(hfs_y1+1:n-hfs_y2,hfs_x1+1:m-hfs_x2);
% 
% end%function
% 
% 
% 
% function x = deconvL2_w(I,filt1,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy)
% 
% if (~exist('max_it','var'))
%    max_it=200;
% end
% 
% [n,m]=size(I);
% 
% hfs1_x1=floor((size(filt1,2)-1)/2);
% hfs1_x2=ceil((size(filt1,2)-1)/2);
% hfs1_y1=floor((size(filt1,1)-1)/2);
% hfs1_y2=ceil((size(filt1,1)-1)/2);
% 
% hfs_x1=hfs1_x1;
% hfs_x2=hfs1_x2;
% hfs_y1=hfs1_y1;
% hfs_y2=hfs1_y2;
% 
% m=m+hfs_x1+hfs_x2;
% n=n+hfs_y1+hfs_y2;
% 
% mask=zeros(n,m);
% mask(hfs_y1+1:n-hfs_y2,hfs_x1+1:m-hfs_x2)=1;
% 
% if (~exist('weight_x','var'))
%   weight_x=ones(n,m-1);
%   weight_y=ones(n-1,m);
%   weight_xx=zeros(n,m-2);
%   weight_yy=zeros(n-2,m);
%   weight_xy=zeros(n-1,m-1);
% end
% 
% tI=I;
% x=tI([ones(1,hfs_y1),1:end,end*ones(1,hfs_y2)],[ones(1,hfs_x1),1:end,end*ones(1,hfs_x2)]);
% 
% b=conv2(x.*mask,filt1,'same');
% 
% dxf=[1 -1];
% dyf=[1;-1];
% dyyf=[-1; 2; -1];
% dxxf=[-1, 2, -1];
% dxyf=[-1 1;1 -1];
% 
% if (max(size(filt1)<25))
%   Ax=conv2(conv2(x,rot90(filt1,2),'same').*mask,  filt1,'same');
% else
%   Ax=fftconv(fftconv(x,rot90(filt1,2),'same').*mask,  filt1,'same');
% end
% 
% Ax=Ax+we*conv2(weight_x.*conv2(x,rot90(dxf,2),'valid'),dxf);
% Ax=Ax+we*conv2(weight_y.*conv2(x,rot90(dyf,2),'valid'),dyf);
% Ax=Ax+we*(conv2(weight_xx.*conv2(x,rot90(dxxf,2),'valid'),dxxf));
% Ax=Ax+we*(conv2(weight_yy.*conv2(x,rot90(dyyf,2),'valid'),dyyf));
% Ax=Ax+we*(conv2(weight_xy.*conv2(x,rot90(dxyf,2),'valid'),dxyf));
% 
% r = b - Ax;
% for iter = 1:max_it  
%      rho = (r(:)'*r(:));
%      if ( iter > 1 ),                       % direction vector
%         beta = rho / rho_1;
%         p = r + beta*p;
%      else
%         p = r;
%      end
%      if (max(size(filt1)<25))
%        Ap=conv2(conv2(p,rot90(filt1,2),'same').*mask,  filt1,'same');
%      else  
%        Ap=fftconv(fftconv(p,rot90(filt1,2),'same').*mask,  filt1,'same');
%      end
% 
%      Ap=Ap+we*conv2(weight_x.*conv2(p,rot90(dxf,2),'valid'),dxf);
%      Ap=Ap+we*conv2(weight_y.*conv2(p,rot90(dyf,2),'valid'),dyf);
%      Ap=Ap+we*(conv2(weight_xx.*conv2(p,rot90(dxxf,2),'valid'),dxxf));
%      Ap=Ap+we*(conv2(weight_yy.*conv2(p,rot90(dyyf,2),'valid'),dyyf));
%      Ap=Ap+we*(conv2(weight_xy.*conv2(p,rot90(dxyf,2),'valid'),dxyf));
% 
% 
%      q = Ap;
%      alpha = rho / (p(:)'*q(:) );
%      x = x + alpha * p;                    % update approximation vector
% 
%      r = r - alpha*q;                      % compute residual
% 
%      rho_1 = rho;
% end
% 
% end%function
