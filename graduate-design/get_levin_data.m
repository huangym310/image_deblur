function [B, L] = get_levin_data(opt)

B = zeros(opt.img_size(1), opt.img_size(2), opt.img_num);

for i = 1:opt.img_num
    img_path = [opt.img_name_set num2str(i) '.mat'];
    load(img_path);
    B(:,:,i) = y.*255;
    L = x;
end


end

