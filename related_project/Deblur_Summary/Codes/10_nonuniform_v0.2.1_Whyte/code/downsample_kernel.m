function kernel_s = downsample_kernel(kernel_s,theta_grid_s,kernel_up,theta_grid_up,scale_ratio_k,theta_step,s)
% 	downsample_kernel   Downsample a blur kernel from one scale to a lower scale
% 		kernel_s = downsample_kernel(kernel_s,theta_grid_s,kernel_up,theta_grid_up,scale_ratio_k,theta_step,s)
%
% 		Inputs:
% 				kernel_s        	kernel at scale s (will be replaced in output)
% 				theta_grid_s    	3 x 1 cell array of angles covered by kernel at scale s
% 				kernel_up       	kernel at higher scale
%				theta_grid_up   	3 x 1 cell array of angles covered by kernel at higher scale
%				scale_ratio_k   	ratio between sizes of subsequent scales
%				theta_step      	theta spacing at finest scale
%				s					current scale s
% 
% 		Outputs:
% 				kernel_s			downsampled kernel 
% 
% 	Author:			Oliver Whyte <oliver.whyte@ens.fr>
% 	Date:			August 2010
% 	Copyright:		2010, Oliver Whyte
% 	Reference:		O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
% 	URL:			http://www.di.ens.fr/~whyte/deblurring/

ksize_s = ones(1,3);
ksize_s(1:ndims(kernel_s)) = size(kernel_s);
nsdims_s = ksize_s > 1;
ksize_up = ones(1,3);
ksize_up(1:ndims(kernel_up)) = size(kernel_up);
nsdims_up = ksize_up > 1;
theta_step_s = theta_step/scale_ratio_k^(s-1);
% If going from 2d to 1d, or 3d to 2d, pad smaller kernel to allow use of
% higher dimension interpolation function
for sdim = find(nsdims_up & ~nsdims_s)
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
k_up = squeeze(kernel_up);
if nnz(nsdims_up) == 3
	k_s_full = interp3(ttx_up,tty_up,ttz_up,k_up,ttx_s,tty_s,ttz_s,'linear',0);
elseif nnz(nsdims_up) == 2
	if ksize_up(1) == 1
		k_s_full = interp2(ttx_up,ttz_up,k_up,ttx_s,ttz_s,'*linear',0);
	elseif ksize_up(2) == 1
		k_s_full = interp2(tty_up,ttz_up,k_up,tty_s,ttz_s,'*linear',0);
	else % ksize_up(3) == 1
		k_s_full = interp2(ttx_up,tty_up,k_up,ttx_s,tty_s,'*linear',0);
	end
else % nnz(nsdims_up) == 1
	if ksize_up(1) > 1
		k_s_full = interp1(tty_up,k_up,tty_s,'linear',0);
	elseif ksize_up(2) > 1
		k_s_full = interp1(ttx_up,k_up,ttx_s,'linear',0);
	else % ksize_up(3) > 1
		k_s_full = interp1(ttz_up,k_up,ttz_s,'linear',0);
	end
end

% If going from 2d to 1d, or 3d to 1d, depad kernel after use of corresponding interpolation function
for sdim = find(nsdims_up & ~nsdims_s)
	% Extract central slice (i.e. slice 2 of 3) along dimension sdim
	switch sdim
	case 1
		k_s_full = k_s_full(2,:,:);
	case 2
		k_s_full = k_s_full(:,2,:);
	case 3
		k_s_full = k_s_full(:,:,2);
	end
end

kernel_s = reshape(k_s_full,size(kernel_s));