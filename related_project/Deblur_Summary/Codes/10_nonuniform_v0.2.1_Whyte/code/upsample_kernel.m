function kernel_up = upsample_kernel(k_s_full,theta_grid_s,kernel_up,theta_grid_up,theta_step_s)
% 	upsample_kernel   Upsample a blur kernel from one scale to a higher scale
% 		kernel_up = upsample_kernel(k_s_full,theta_grid_s,kernel_up,theta_grid_up,theta_step_s)
%
% 		Inputs:
% 				k_s_full        	kernel at scale s
% 				theta_grid_s    	3 x 1 cell array of angles covered by kernel at scale s
% 				kernel_up       	kernel at higher scale
%				theta_grid_up   	3 x 1 cell array of angles covered by kernel at higher scale
%				theta_step_s      	theta spacing at current scale s
% 
% 		Outputs:
% 				kernel_up			upsampled kernel 
% 
% 	Author:			Oliver Whyte <oliver.whyte@ens.fr>
% 	Date:			August 2010
%	Copyright:		2010, Oliver Whyte
% 	Reference:		O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
% 	URL:			http://www.di.ens.fr/~whyte/deblurring/

for i=1:3
    theta_grid_s{i} = reshape(theta_grid_s{i},size(k_s_full));
    theta_grid_up{i} = reshape(theta_grid_up{i},size(kernel_up));
end

ksize_s = ones(1,3);
ksize_s(1:ndims(k_s_full)) = size(k_s_full);
nsdims_s = ksize_s > 1;
ksize_up = ones(1,3);
ksize_up(1:ndims(kernel_up)) = size(kernel_up);
nsdims_up = ksize_up > 1;
% If going from 1d to 2d, or 2d to 3d, pad kernel to allow use of corresponding interpolation function
for sdim = find(nsdims_up & ~nsdims_s)
    % Pad kernel along singular dimension with a layer of zeros either side
    k_s_full = cat(sdim,zeros(size(k_s_full)),k_s_full,zeros(size(k_s_full)));
    % Extend theta arrays:
    % For singular dimension, extrapolate values
    theta_grid_s{sdim} = cat(sdim,theta_grid_s{sdim}-theta_step_s(sdim),theta_grid_s{sdim},theta_grid_s{sdim}+theta_step_s(sdim));
    % For other dimensions, replicate values
    repmatdims = [1, 1, 1];      repmatdims(sdim) = 3;
    for nsdim = find([1, 2, 3] ~= sdim)
        theta_grid_s{nsdim} = repmat(theta_grid_s{nsdim},repmatdims);
    end
end
tty_s = squeeze(theta_grid_s{1});
ttx_s = squeeze(theta_grid_s{2});
ttz_s = squeeze(theta_grid_s{3});
tty_up = squeeze(theta_grid_up{1});
ttx_up = squeeze(theta_grid_up{2});
ttz_up = squeeze(theta_grid_up{3});
k_s_full = squeeze(k_s_full);
if nnz(nsdims_up) == 3
    k_up = interp3(ttx_s,tty_s,ttz_s,k_s_full,ttx_up,tty_up,ttz_up,'linear',0);
elseif nnz(nsdims_up) == 2
    if ksize_up(1) == 1
        k_up = interp2(ttx_s,ttz_s,k_s_full,ttx_up,ttz_up,'*linear',0);
    elseif ksize_up(2) == 1
        k_up = interp2(tty_s,ttz_s,k_s_full,tty_up,ttz_up,'*linear',0);
    else % ksize_up(3) == 1
        k_up = interp2(ttx_s,tty_s,k_s_full,ttx_up,tty_up,'*linear',0);
    end
else % nnz(nsdims_up) == 1
    if ksize_up(1) > 1
        k_up = interp1(tty_s,k_s_full,tty_up,'linear',0);
    elseif ksize_up(2) > 1
        k_up = interp1(ttx_s,k_s_full,ttx_up,'linear',0);
    else % ksize_up(3) > 1
        k_up = interp1(ttz_s,k_s_full,ttz_up,'linear',0);
    end
end
kernel_up = reshape(k_up,size(kernel_up));
