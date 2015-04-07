function K = initial_kernel(opt)

K = zeros(opt.ker_size(1), opt.ker_size(2), opt.img_num);
% K(opt.ker_size(1),opt.ker_size(2),:) = 0.6;
% K(opt.ker_size(1) - 1, opt.ker_size(2),:) = 0.2;
% K(opt.ker_size(1),opt.ker_size(2) - 1,:) = 0.2;
K = 1/(opt.ker_size(1) * opt.ker_size(2)) * ones(opt.ker_size(1), opt.ker_size(2), opt.img_num);
% K = zeros(opt.ker_size(1),opt.ker_size(2),opt.img_num);
% K(8,3:12,1) = 0.05;
% K(9,3:12,1) = 0.05;
% K(3:12,8,2) = 0.05;
% K(3:12,9,2) = 0.05;

end

