addpath Levin_data
addpath Levin

num_image = 4;
ker_size = 19;
for i = 1:4
    K = zeros(ker_size, ker_size, num_image);
    for j = 1:num_image
        load(['im0' num2str(i) '_ker0' num2str(j) '.mat']);
        if j == 1
            img_size = size(x);
            B = zeros(img_size(1) - ker_size + 1, img_size(2) - ker_size + 1, num_image);
            %B = zeros(size(L,1), size(L,2), num_image);
            L = x;
        end
        K(:,:,j) = imresize(f, [ker_size ker_size]);
        K(:,:,j) = K(:,:,j)./sum(sum(K(:,:,j)));
        B(:,:,j) = conv2(L, K(:,:,j), 'valid');
    end
    save(['Levin/im0' num2str(i) '.mat'],'L','K','B');
end