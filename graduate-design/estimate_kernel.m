function K = estimate_kernel(B, opt)
%% blind kernel estimation from multiple image.
K = initial_kernel(opt);

Dkx = zeros(opt.ker_size(1)*(opt.ker_size(2) - 1)*opt.img_num, opt.ker_size(1)*opt.ker_size(2)*opt.img_num);
Dky = zeros(opt.ker_size(1)*(opt.ker_size(2) - 1)*opt.img_num, opt.ker_size(1)*opt.ker_size(2)*opt.img_num);
base_row = 0;  base_col = 0;
for i = 1:opt.img_num
    for j = 1:opt.ker_size(1)*(opt.ker_size(2) - 1)
        Dkx(base_row + j, base_col + j)  = -1;
        Dkx(base_row + j, base_col + opt.ker_size(1) + j) = 1;
    end
    base_row = base_row + opt.ker_size(1) * (opt.ker_size(2) - 1);
    base_col = base_col + opt.ker_size(1) * opt.ker_size(2);
end
base_row = 0;  base_col = 0;
for i = 1:opt.img_num
    for j = 1:opt.ker_size(1)
        for k = 1:opt.ker_size(2) - 1
            Dky(base_row + k, base_col + k) = -1;
            Dky(base_row + k, base_col + k + 1) = 1;
        end
        base_row = base_row + opt.ker_size(2) - 1;
        base_col = base_col + opt.ker_size(2);
    end
end

dx = [-1,1;0,0];
dy = [-1,0;1,0];
LoG = [0,0,1,0,0;0,1,2,1,0;1,2,-16,2,1;0,1,2,1,0;0,0,1,0,0];
% Bx = B;
% By = B;
Bx = zeros(opt.img_size(1) - 4, opt.img_size(2) - 4, opt.img_num);
for i = 1:opt.img_num
    Bx(:,:,i) = conv2(B(:,:,i), LoG,'valid');
end
% Bx = B;
% By = B;
% Bx = zeros(opt.img_size(1) - 1, opt.img_size(2) - 1, opt.img_num);
% By = zeros(opt.img_size(1) - 1, opt.img_size(2) - 1, opt.img_num);
% for i = 1:opt.img_num
%     Bx(:,:,i) = conv2(B(:,:,i),dx,'valid');
%     By(:,:,i) = conv2(B(:,:,i),dy,'valid');
% end

img_size = size(Bx(:,:,1));
% Bx_pad = wrap_boundary(Bx,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);
% By_pad = wrap_boundary(By,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);

% the step below need to be modified to get avoid of the explicit storage
% of G if the memory space is limited.
% Gx = cell(opt.img_num,1);
% Gy = cell(opt.img_num,1);
% for i = 1:opt.img_num
%     Gx{i} = Get_Toeplitz(Bx(:,:,i), img_size, opt.ker_size);
%     Gy{i} = Get_Toeplitz(By(:,:,i), img_size, opt.ker_size);
% end
Gx = cell(opt.img_num,1);
for i = 1:opt.img_num
    Gx{i} = Get_Toeplitz(Bx(:,:,i), img_size, opt.ker_size);
end

block_size = size(Gx{1},2);
R = zeros(block_size * opt.img_num, block_size * opt.img_num);
tr_mat = zeros(block_size, block_size);
for i = 1:opt.img_num
    for j = 1:opt.img_num
        R((i-1)*block_size + 1: i*block_size, (j-1)*block_size + 1: j*block_size) = - Gx{j}' * Gx{i};
        if i == j
            tr_mat = tr_mat + Gx{i}'* Gx{i};
        end
%         R((i-1)*block_size + 1: i*block_size, (j-1)*block_size + 1: j*block_size) = - Gx{j}' * Gx{i} - Gy{j}' * Gy{i};
%         if i == j
%             tr_mat = tr_mat + Gx{i}'* Gx{i} + Gy{i}'* Gy{i};
%         end
    end
end

for i = 1:opt.img_num
    R((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) = ...
        R((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) + tr_mat;
end
clear GX 

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

W = zeros(length(k),length(k));
% fixed penalty coefficiency. May need to modify later.
iter = 0;
while ~converge
    
   % update sparse Weight matrix W;
   for i = 1:length(k)
       W(i,i) = 1/(abs(k(i)^(2 - opt.p)) + 0.000001);
   end
   
   % update smooth weight matrix Wx and Wy
   
   M = 2 * R + opt.sparse.* W + opt.smooth.*(Dkx' * Dkx + Dky' * Dky) + mu * eye(opt.img_num * opt.ker_size(1) * opt.ker_size(2));
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
   if max(norm(x - x_prev), norm(k - k_prev)) < epsilon
        mu = min(mu * beta, mu_max);
   end
   
   % convergence condition
   cri = max([norm(k - x), norm(x - x_prev), norm(k - k_prev)]);
   disp(cri);
   if  iter > opt.max_iter
       converge = 1;
   end
   
   x_prev = x;
   k_prev = k;
   iter = iter + 1;
   
end

estimation_error = 2 * k' * R * k;
sparse_penalty = k'* W * k;
smooth_penalty = k' * (Dkx' * Dkx + Dky' * Dky) * k;

K = reshape(x,[opt.ker_size(1),opt.ker_size(2), opt.img_num]);

end

