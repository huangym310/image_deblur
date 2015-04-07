function K = estimate_kernel_v6(B, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     blind kernel estimation from multiple image.
%%
%% use kernel prior proposed by Liu. et.al in TIP2014
%%
%%           Blind Image Deblurring Using
%%   Spectral Properties of Convolution Operators
%%
%%    + the pyramid of image comparing to v5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda1 = 0;
lambda2 = 1;

K1 = zeros(15,15);
K2 = zeros(15,15);
K1(8,3:12) = 0.1;
K2(3:12,8) = 0.1;

ker_scale = [1;3;5;7;9;11;13;15;17;19;23];
K = zeros(ker_scale(1),ker_scale(1),opt.img_num);

LoG = [0,0,1,0,0;0,1,2,1,0;1,2,-16,2,1;0,1,2,1,0;0,0,1,0,0];
for i = 1:opt.img_num
    Bx(:,:,i) = conv2(B(:,:,i), LoG, 'same');
end

for s = 2:length(ker_scale)-1
    
    % the step below need to be modified to get avoid of the explicit storage
    % of G if the memory space is limited.
    G = cell(opt.img_num,1);
    for i = 1:opt.img_num
        G{i} = conv2mtx(Bx(:,:,i), ker_scale(s), ker_scale(s),'same');
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
        G{i} = conv2mtx(Bx(:,:,i), ker_scale(s+1), ker_scale(s+1), 'same');
        M{i} = G{i}' * G{i};
        [E{i},V{i}] = eig(M{i});
        H{i} = zeros(ker_scale(s) * ker_scale(s), ker_scale(s) * ker_scale(s));
        for j = 1:size(E{i},1)
            kappa = reshape(E{i}(:,j), [ker_scale(s+1) ker_scale(s+1)]);
            % kappa = reshape(E{i}(:,j), opt.ker_size);
            % Toeplitz_kapp = conv2mtx(kappa, 1, 1, 'valid')';
            Toeplitz_kapp = conv2mtx(kappa, ker_scale(s), ker_scale(s), 'same');
            H{i} = H{i} + (Toeplitz_kapp'*Toeplitz_kapp)/V{i}(j,j);
        end
        R2((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) = ...
            R2((i-1)*block_size + 1: i*block_size, (i-1)*block_size + 1:i*block_size) + H{i};
    end
    
    R = lambda1 * R1 + lambda2 * R2;
    K = imresize(K,[ker_scale(s), ker_scale(s)]);
    k = K(:);
    e = ones(opt.img_num,1);
    C = zeros(opt.img_num, opt.img_num* ker_scale(s)* ker_scale(s));
    for i = 1:opt.img_num
        C(i, (i-1)* ker_scale(s)* ker_scale(s) + 1 : i * ker_scale(s)* ker_scale(s)) = ones(1, ker_scale(s) * ker_scale(s));
    end
    mu = opt.mu_init;
    mu_max = opt.mu_max;
    beta = opt.beta;
    epsilon = opt.epsilon;
    
    x = zeros(opt.img_num * ker_scale(s) * ker_scale(s), 1);
    y = zeros(opt.img_num * ker_scale(s) * ker_scale(s), 1);
    converge = 0;
    k_prev = zeros(length(k),1);
    x_prev = zeros(length(x),1);
    
    % fixed penalty coefficiency. May need to modify later.
    iter = 0;
    while ~converge
        
        M = 2 * R + mu * eye(opt.img_num * ker_scale(s) * ker_scale(s));
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
        
        if mod(iter,50) == 0
            K = reshape(x, [ker_scale(s), ker_scale(s), opt.img_num]);
            figure(1);      subplot(2,2,1);    imshow(K1,[]);
            figure(1);      subplot(2,2,2);    imshow(K2,[]);
            figure(1);      subplot(2,2,3);    imshow(K(:,:,1),[]);
            figure(1);      subplot(2,2,4);    imshow(K(:,:,2),[]);
            pause(0.05);
        end
        
        x_prev = x;
        k_prev = k;
        iter = iter + 1;
    end
    K = reshape(x, [ker_scale(s), ker_scale(s), opt.img_num]);
    
end


end

