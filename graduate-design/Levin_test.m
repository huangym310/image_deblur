clear;
addpath Levin
load('im04.mat');

opt.img_size = size(B(:,:,1));
opt.ker_size = size(K(:,:,1));
opt.img_num = 4;

opt.beta = 1.05;
opt.mu_init = 0.00001;
opt.mu_max = 10000;
opt.lambda = 0.1;
opt.max_iter = 200;
opt.epsilon = 1e-8;
opt.sparse = 1;
opt.smooth = 0.15;
opt.p = 0.4;

Ks = estimate_kernel_v8(B,K,opt);
S = zeros(size(B));
for i = 1:opt.img_num
    S(:,:,i)   = nonblind_deconv(B(:,:,i),K(:,:,i),1e-5);
    Ss(:,:,i)  = nonblind_deconv(B(:,:,i),Ks(:,:,i),1e-5);
    figure(2);    
    subplot(4,opt.img_num, i);                   imshow(B(:,:,i),[min(min(L)), max(max(L))]);
    subplot(4,opt.img_num, opt.img_num + i);     imshow(S(:,:,i),[min(min(L)), max(max(L))]);   
    subplot(4,opt.img_num, 2*opt.img_num + i);   imshow(Ss(:,:,i),[min(min(L)), max(max(L))]);
    subplot(4,opt.img_num, 3*opt.img_num + i);   imshow(L,[min(min(L)), max(max(L))]);
end
% 
S_mul  = nonblind_deconv(B,K,1e-5);
Ss_mul = nonblind_deconv(B,Ks,1e-5);
% 
figure(3);
subplot(1,4,1);  imshow(S_mul,[min(min(L)),  max(max(L))]);      title('deblur with ground truth kernel');
subplot(1,4,2);  imshow(Ss_mul,[min(min(L)), max(max(L))]);      title('deblur with estimated kernel');
subplot(1,4,4);  imshow(L,[min(min(L)),      max(max(L))]);      title('Ground truth image');

D = zeros(size(Ss,3), size(Ss,1) * size(Ss,2));
for i = 1:opt.img_num
    D(i,:) = vec(Ss(:,:,i));
end

[A,E] = rpca_adm(D,0.02);
A = reshape(A(1,:),[size(B,1), size(B,2)]);
figure(3)
subplot(1,4,3);   imshow(A,[min(min(L)), max(max(L))]);         title('deblur with RPCA');
