function K = estimate_kernel_v5(B, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     blind kernel estimation from multiple image.
%%
%% use kernel prior proposed by Liu. et.al in TIP2014
%%
%%           Blind Image Deblurring Using 
%%   Spectral Properties of Convolution Operators
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda1 = 0;
lambda2 = 1;

%B_pad = wrap_boundary(B,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]); 
K = initial_kernel(opt);
LoG = [0,0,1,0,0;0,1,2,1,0;1,2,-16,2,1;0,1,2,1,0;0,0,1,0,0];
% Bx = B;
% By = B;
% Bx = zeros(opt.img_size(1), opt.img_size(2), opt.img_num);
for i = 1:opt.img_num
%    Bx(:,:,i) = B(:,:,i);
    Bx(:,:,i) = conv2(B(:,:,i), LoG, 'valid');
end

% Bx_pad = wrap_boundary(Bx,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);
% By_pad = wrap_boundary(By,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);

% the step below need to be modified to get avoid of the explicit storage
% of G if the memory space is limited.
G = cell(opt.img_num,1);

for i = 1:opt.img_num
    G{i} = conv2mtx(Bx(:,:,i), opt.ker_size(1), opt.ker_size(2),'valid');
end

block_size = size(G{1},2);
R = zeros(block_size * opt.img_num, block_size * opt.img_num);
R1 = zeros(block_size * opt.img_num, block_size * opt.img_num);
R2 = zeros(block_size * opt.img_num, block_size * opt.img_num);
tr_mat = zeros(block_size, block_size);
for i = 1:opt.img_num
    for j = 1:opt.img_num
        R1((i-1)*block_size + 1: i*block_size, (j-1)*block_size + 1: j*block_size) = - G{j}' * G{i};
        if i == j
            tr_mat = tr_mat + G{i}'* G{i};
        end
    end
end

for i = 1:opt.img_num
    R1((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) = ...
        R1((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) + tr_mat;
end

clear G

%% calculating the kernel prior matrix H
M = cell(opt.img_num,1);
E = cell(opt.img_num,1);
V = cell(opt.img_num,1);
H = cell(opt.img_num,1);
for i = 1:opt.img_num
     %G{i} = conv2mtx(Bx(:,:,i), opt.ker_size(1), opt.ker_size(2),'valid');
     G{i} = conv2mtx(Bx(:,:,i), opt.smp_size(1), opt.smp_size(2),'valid');
     M{i} = G{i}' * G{i};
     [E{i},V{i}] = eig(M{i});
     H{i} = zeros(opt.ker_size(1) * opt.ker_size(2), opt.ker_size(1) * opt.ker_size(2));
     for j = 1:size(E{i},1)
         kappa = reshape(E{i}(:,j), opt.smp_size);
         % kappa = reshape(E{i}(:,j), opt.ker_size);
         % Toeplitz_kapp = conv2mtx(kappa, 1, 1, 'valid')';
         Toeplitz_kapp = conv2mtx(kappa, opt.ker_size(1), opt.ker_size(2), 'valid');
         H{i} = H{i} + (Toeplitz_kapp'*Toeplitz_kapp)/V{i}(j,j); 
     end
     R2((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) = ...
        R2((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) + H{i};
end

R = lambda1 * R1 + lambda2 * R2;
   
k = K(:);
e = ones(opt.img_num,1);
C = zeros(opt.img_num, opt.img_num* opt.ker_size(1)* opt.ker_size(2));
for i = 1:opt.img_num
    C(i, (i-1)* opt.ker_size(1)* opt.ker_size(2) + 1 : i * opt.ker_size(1)* opt.ker_size(2)) = ones(1, opt.ker_size(1) * opt.ker_size(2));  
end
mu = opt.mu_init;
mu_max = opt.mu_max;
beta = opt.beta;
epsilon = opt.epsilon;

x = zeros(opt.img_num * opt.ker_size(1) * opt.ker_size(2),1);
y = zeros(opt.img_num * opt.ker_size(1) * opt.ker_size(2),1);
converge = 0;
k_prev = zeros(length(k),1);
x_prev = zeros(length(x),1);

% fixed penalty coefficiency. May need to modify later.
iter = 0;
while ~converge
   
   M = 2 * R + mu * eye(opt.img_num * opt.ker_size(1) * opt.ker_size(2));
   tmp = inv(M);

   % update the kernel k
   lambda = (C * tmp * C')\(C * tmp * (mu * x - y) - e); 
   k = tmp * (mu * x - y) - tmp * C' * lambda;
   
   % update x to ensure the nonnegativeness of k.
   for i = 1:length(x)
       x(i) = max(0, k(i)+ y(i)/mu);
   end
   
   % update dual variable.
   y = y + mu*(k - x);
   
   % update mu
   if max(norm(x - x_prev), norm(k - k_prev)) < opt.epsilon
        mu = min(mu * beta, mu_max);
   end
   
   % convergence condition
   cri = max([norm(k - x), norm(x - x_prev), norm(k - k_prev)]);
   disp(cri);
   %if  cri < opt.epsilon
   if iter == 200
       converge = 1;
   end
   
   x_prev = x;
   k_prev = k;
   iter = iter + 1;
   
end

K = reshape(x,[opt.ker_size(1),opt.ker_size(2), opt.img_num]);

end

