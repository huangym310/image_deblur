function B = prepare_data(opt)

B = zeros(opt.img_size(1), opt.img_size(2), opt.img_num);

for i = 1:opt.img_num
    img_path = [opt.img_name_set num2str(i) opt.img_format];
    I = imread(img_path);
    if size(I,3) ~= 1
        I = rgb2gray(I);
    end
    I = im2double(I);
    I = imresize(I, opt.img_size);
    B(:,:,i) = I;
end

