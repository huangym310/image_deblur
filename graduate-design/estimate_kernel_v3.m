function K = estimate_kernel_v3(B, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% blind kernel estimation from multiple image.
%%
%% + LoG filter for preprocessing comparing to v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% configure the para in AM algorithm
para.mu1 = 1000;
para.mu2 = 1;
para.mu3 = 1;
para.mu4 = 1000;
para.lam1 = 0.00015;
para.lam2 = 0.00003;

%% compute the matrix Dkx and Dky 
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

%% compute matrix R so that ||B2 * k1 - B1 * k2|| = k' * R * k
LoG = [0,0,1,0,0;0,1,2,1,0;1,2,-16,2,1;0,1,2,1,0;0,0,1,0,0];
% Bx = B;
% By = B;
Bx = zeros(opt.img_size(1) - 4, opt.img_size(2) - 4, opt.img_num);
for i = 1:opt.img_num
    Bx(:,:,i) = conv2(B(:,:,i), LoG,'valid');
end

img_size = size(Bx(:,:,1));
% Bx_pad = wrap_boundary(Bx,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);
% By_pad = wrap_boundary(By,[img_size(1)+opt.ker_size(1) - 1  img_size(2)+opt.ker_size(2) - 1]);

% the step below need to be modified to get avoid of the explicit storage
% of G if the memory space is limited.
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
    end
end
for i = 1:opt.img_num
    R((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) = ...
        R((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) + tr_mat;
end
clear GX 

%% equality constraint on k so that |k1| = 1 && |k2| = 1  <=>  Ck = e
e = ones(opt.img_num,1);
C = zeros(opt.img_num, opt.img_num* opt.ker_size(1)* opt.ker_size(2));
for i = 1:opt.img_num
    C(i, (i-1)* opt.ker_size(1)* opt.ker_size(2) + 1 : i * opt.ker_size(1)* opt.ker_size(2)) = ones(1, opt.ker_size(1) * opt.ker_size(2));  
end

K = initial_kernel(opt);
k = K(:);
n = length(k);
w = k;
u = k;
gx = Dkx * k;
gy = Dky * k;
prev_k = k;
prev_w = w;
prev_u = u;
prev_gx = gx;
prev_gy = gy;

A = R + para.mu1 * C' * C + (para.mu2 + para.mu4)* eye(n,n) + para.mu3 * (Dkx' * Dkx + Dky' * Dky);
inv_A = inv(A);
converge = 0;
iter = 0;
inner_iter = 0;
thre = 0.001;
while ~converge
    b =  2*(para.mu1 * e' * C + para.mu2 * u' + para.mu3 * (gx'* Dkx + gy' *Dky) + para.mu4 * w');
    k = (1/2) * inv_A * b';
    
    for i = 1:length(k)
       if k(i)^2 >= para.lam1/para.mu2
          u(i) = k(i); 
       else 
          u(i) = 0;
       end
    end
    
    w = max(k,0);
    
    kx = Dkx * k;
    ky = Dky * k;
    for i = 1:length(kx)
        if kx(i)^2 >= para.lam2 / para.mu3
            gx(i) = kx(i);
        else
            gx(i) = 0;
        end
    end
    for i = 1:length(ky)
        if ky(i)^2 >= para.lam2 / para.mu3
            gy(i) = ky(i);
        else
            gy(i) = 0;
        end
    end
    
    varia = max(abs([norm(k-prev_k,inf), norm(u-prev_u,inf), norm(w-prev_w,inf), ...
                    norm(gx-prev_gx,inf), norm(gy-prev_gy,inf)]));
    resi = max(abs([norm(C*k -e, inf), norm(k - u, inf), norm(k - w, inf),...
                    norm(kx - gx, inf), norm(ky - gy, inf)]));
    cri = max(varia,resi);
    disp(cri);
    if cri < thre || inner_iter == 50
        inner_iter = 0;
        iter = iter + 1;
        para.mu1 = para.mu1 * 1.01;
        para.mu2 = para.mu2 * 1.01;
        para.mu3 = para.mu3 * 1.01;
        para.mu4 = para.mu4 * 1.01;
        A = R + para.mu1 * C' * C + (para.mu2 + para.mu4)* eye(n,n) + para.mu3 * (Dkx' * Dkx + Dky' * Dky);
        inv_A = inv(A);
    else 
        inner_iter = inner_iter + 1;
    end
    
    if iter == 500
        break;
    end

    prev_k = k; 
    prev_u = u;
    prev_w = w;
    prev_gx = gx;
    prev_gy = gy;
end

K = reshape(k,[opt.ker_size(1),opt.ker_size(2), opt.img_num]);