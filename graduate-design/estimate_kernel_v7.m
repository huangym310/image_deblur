function K = estimate_kernel_v7(B,K_grd,opt)
%% blind kernel estimation from multiple image.

ker_scale = [1;3;5;7;9;11;13;15;17;19];
K = ones(ker_scale(1), ker_scale(1), opt.img_num)./(ker_scale(1) * ker_scale(1));

for s = 2:length(ker_scale)
    opt.sparse = opt.sparse * 0.9;
    opt.smooth = opt.smooth * 0.89;
    Bx = imresize(B, ker_scale(s)/ker_scale(end));
    LoG = [0,0,1,0,0;0,1,2,1,0;1,2,-16,2,1;0,1,2,1,0;0,0,1,0,0];
    for i = 1:opt.img_num
        Bx(:,:,i) = conv2(Bx(:,:,i), LoG, 'same');
    end
    img_size = size(Bx(:,:,1));
    
    Gx = cell(opt.img_num,1);
    for i = 1:opt.img_num
        Gx{i} = Get_Toeplitz(Bx(:,:,i), img_size, [ker_scale(s) ker_scale(s)]);
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
    
    Dkx = zeros(ker_scale(s) * (ker_scale(s) - 1) * opt.img_num, ker_scale(s) * ker_scale(s) * opt.img_num);
    Dky = zeros(ker_scale(s) * (ker_scale(s) - 1) * opt.img_num, ker_scale(s) * ker_scale(s) * opt.img_num);
    base_row = 0;  base_col = 0;
    for i = 1:opt.img_num
        for j = 1 : ker_scale(s) * (ker_scale(s) - 1)
            Dkx(base_row + j, base_col + j)  = -1;
            Dkx(base_row + j, base_col + ker_scale(s) + j) = 1;
        end
        base_row = base_row + ker_scale(s) * (ker_scale(s) - 1);
        base_col = base_col + ker_scale(s) * ker_scale(s);
    end
    base_row = 0;  base_col = 0;
    for i = 1:opt.img_num
        for j = 1:ker_scale(s)
            for k = 1:ker_scale(s) - 1
                Dky(base_row + k, base_col + k) = -1;
                Dky(base_row + k, base_col + k + 1) = 1;
            end
            base_row = base_row + ker_scale(s) - 1;
            base_col = base_col + ker_scale(s);
        end
    end
    K = imresize(K, [ker_scale(s), ker_scale(s)]);
    for i = 1:opt.img_num
        K(:,:,i) = K(:,:,i)./sum(sum(K(:,:,i)));
    end
    k = K(:);
    e = ones(opt.img_num,1);
    C = zeros(opt.img_num, opt.img_num* ker_scale(s)* ker_scale(s));
    for i = 1:opt.img_num
        C(i, (i-1)* ker_scale(s) * ker_scale(s) + 1 : i * ker_scale(s)* ker_scale(s)) = ones(1, ker_scale(s) * ker_scale(s));
    end
    mu = opt.mu_init;
    mu_max = opt.mu_max;
    beta = opt.beta;
    epsilon = opt.epsilon;
    
    x = zeros(opt.img_num * ker_scale(s) * ker_scale(s),1);
    y = zeros(opt.img_num * ker_scale(s) * ker_scale(s),1);
    converge = 0;
    k_prev = zeros(length(k),1);
    x_prev = zeros(length(x),1);
    
    W = zeros(length(k),length(k));
    Wx = zeros(ker_scale(s) * (ker_scale(s)-1) * opt.img_num, ker_scale(s) * (ker_scale(s)-1) * opt.img_num);
    Wy = zeros(ker_scale(s) * (ker_scale(s)-1) * opt.img_num, ker_scale(s) * (ker_scale(s)-1) * opt.img_num);
    % fixed penalty coefficiency. May need to modify later.
    iter = 0;
    D = Dkx' * Dkx + Dky' * Dky;
    I = eye(opt.img_num * ker_scale(s)* ker_scale(s));
    while ~converge
        
        % update sparse Weight matrix W;
        for i = 1:length(k)
            W(i,i) = 1/(abs(k(i)^(2 - opt.p)) + 0.00008);
        end
        
        % update smooth weight matrix Wx and Wy
        kx = Dkx * k;
        ky = Dky * k;
        for i = 1:length(kx)
            Wx(i,i) = 1/(abs(kx(i)^(2 - opt.p)) + 0.00008);
            Wy(i,i) = 1/(abs(ky(i)^(2 - opt.p)) + 0.00008);
        end
%             
%         %M = 2 * R + opt.sparse.* W + opt.smooth.*D + mu * I;
        M = 2 * R + opt.sparse.* W + opt.smooth.*(Dkx' * Wx * Dkx + Dky' * Wy * Dky) + mu * I;
        lu = zeros(length(k),1);
        k = quadprog(M,[],[],[],C,e,lu,[]);
%         tmp = inv(M);
%         
%         % update the kernel k
%         lambda = (C * tmp * C')\(C * tmp * (mu * x - y) - e);
%         k = tmp * (mu * x - y) - tmp * C' * lambda;
%         
%         % update x to ensure the nonnegativeness of k.
%         for i = 1:length(x)
%             x(i) = max(0, k(i)+ y(i)/mu);
%         end
%         
%         % update dual variable.
%         y = y + mu*(k - x);
%         
%         % update mu
%         if max(norm(x - x_prev), norm(k - k_prev)) < epsilon
%             mu = min(mu * beta, mu_max);
%         end
%         
        % convergence condition
%         cri = min(max([norm(x - x_prev), norm(k - k_prev)]), norm(k - x));
        cri = norm(k - k_prev, 'inf');
        disp(cri);
        if  iter == 500 || cri < opt.epsilon
            converge = 1;
        end
        
        if mod(iter,100) == 0
            K = reshape(x, [ker_scale(s), ker_scale(s), opt.img_num]);
            for i = 1:opt.img_num
                figure(1);      
                subplot(2,opt.img_num, i);                  imshow(K_grd(:,:,i),[]);
                subplot(2,opt.img_num, opt.img_num + i);    imshow(K(:,:,i),[]);
            end
            pause(0.05);
        end
%         x_prev = x;
        k_prev = k;
        iter = iter + 1;
    end
    
    K = reshape(x,[ker_scale(s),ker_scale(s), opt.img_num]);
    if s == length(ker_scale)
%         test_data_penalty = 2 * k' * R * k;
%         test_sparse_penalty = k'* W * k;
%         test_smooth_penalty = k' * (Dkx' * Dkx + Dky' * Dky) * k;
%         k_grd = K_grd(:);
%         grd_data_penalty = 2 * k_grd' * R * k_grd;
%         grd_sparse_penalty = k_grd'* W * k_grd;
%         grd_smooth_penalty = k_grd' * (Dkx' * Dkx + Dky' * Dky) * k_grd;
        
        K = imresize(K, [opt.ker_size(1), opt.ker_size(2)]);
        for i = 1:opt.img_num
            K(:,:,i) = K(:,:,i)./sum(sum(K(:,:,i)));
        end
    end
end

end

