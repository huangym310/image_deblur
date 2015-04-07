function [pyr_kernel,pyr_tt,pyr_tgs] = make_kernel_pyramid(theta_x_lims,theta_y_lims,theta_z_lims,tgs,scale_ratio_k,max_levels,non_uniform,init_kernel,chain_downsampling,init_tt)
% 	make_kernel_pyramid		Make a scale pyramid for a blur kernel
%		[pyr_kernel,pyr_tt,pyr_tgs] = make_kernel_pyramid(theta_x_lims,theta_y_lims,theta_z_lims,tgs,scale_ratio_k,max_levels,non_uniform,init_kernel,chain_downsampling,init_tt)
%
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/


% sequentially downsample from one scale to the next? otherwise make each
% scale by interpolating the finest one
if nargin < 8, init_kernel = []; end
if nargin < 9, chain_downsampling = false; end
if nargin < 10, init_tt = []; end

% Theta grid spacing
if numel(tgs)==1
    tgs = [tgs, tgs, tgs];
end

%     if non_uniform
% Angles for CSF, not necessarily centred on 0
khh = ceil((theta_y_lims(2)-theta_y_lims(1))/2/tgs(1)-0.5); % half height
khw = ceil((theta_x_lims(2)-theta_x_lims(1))/2/tgs(2)-0.5); % half width
khd = ceil((theta_z_lims(2)-theta_z_lims(1))/2/tgs(3)-0.5); % half depth
ttc = [(theta_y_lims(2)+theta_y_lims(1))/2; ... central theta_y
	   (theta_x_lims(2)+theta_x_lims(1))/2; ... central theta_x
	   (theta_z_lims(2)+theta_z_lims(1))/2]; %  central theta_z
if ~non_uniform
    if ~isequal(ttc,round(ttc))
        error('kernel centre must be a whole number for uniform blur');
    end
end
%     else
% % Offsets for uniform PSF, centred on 0
%         khh = ceil(max(theta_y_lims)/theta_grid_size-0.5); % half height
%         khw = ceil(max(theta_x_lims)/theta_grid_size-0.5); % half width
%         khd = ceil(max(theta_z_lims)/theta_grid_size-0.5); % half depth
%         theta_grid_size = 1; % pixel grid
%         ttc = [0;0;0];
%     end

pyr_kernel = cell(1);
pyr_tt = cell(1,3);
pyr_tgs = cell(1);
s = 1;
while s <= max_levels
    scale_factor = scale_ratio_k^(s-1);
	if non_uniform
	    tgs_s = tgs/scale_factor;
	else
		tgs_s = tgs;
	end
    pyr_tgs{s} = tgs_s;
    % Find size of smallest array which covers same range of angles (or translations) as full-res kernel
    khh_s = ceil(scale_factor*(khh+0.5)-0.5);
    khw_s = ceil(scale_factor*(khw+0.5)-0.5);
    khd_s = ceil(scale_factor*(khd+0.5)-0.5);
    if max([khh_s,khw_s,khd_s]) < 2; max_levels = s; end
    [ttx_s,tty_s,ttz_s] = meshgrid((-khw_s:khw_s)*tgs_s(2)+ttc(2),(-khh_s:khh_s)*tgs_s(1)+ttc(1),(-khd_s:khd_s)*tgs_s(3)+ttc(3));
    pyr_tt{s,1} = tty_s;
    pyr_tt{s,2} = ttx_s;
    pyr_tt{s,3} = ttz_s;
    pyr_kernel{s} = zeros(size(tty_s));
    if ~isempty(init_kernel)
        if s==1
            % If initial kernel specified as string, make the actual kernel
            if ischar(init_kernel)
                if strcmp(init_kernel,'delta')
                    init_kernel =  delta_kernel_rot(2*[khh_s,khw_s,khd_s]+1);
                else
                    error('Unknown string argument for kernel initialisation');
				end
			end
			if isempty(init_tt)
				pyr_kernel{s} = init_kernel;
			else
				% Place a Gaussian at each sample point on the rotation curve, and sum Gaussians at each element in CSF
				% pyr_kernel{s} = griddata3(init_tt(1,:)',init_tt(2,:)',init_tt(3,:)',init_kernel(:),ttx_s,tty_s,ttz_s);
				blur_kernel = zeros(size(ttx_s));
				for k=1:size(init_tt,2)
					disps = [tty_s(:)' - init_tt(2,k); ...
					ttx_s(:)' - init_tt(1,k); ...
					ttz_s(:)' - init_tt(3,k)];
					dispsT = inv(diag((tgs/3).^2))*disps;
					dists_sqr = sum(disps.*dispsT,1);
					blur_kernel = blur_kernel + reshape(exp(-dists_sqr/2),size(ttx_s));
				end	
				pyr_kernel{s} = blur_kernel;
			end
        else
			if non_uniform
				if chain_downsampling
					pyr_kernel{s} = downsample_kernel(pyr_kernel{s},pyr_tt(s,:),pyr_kernel{s-1},pyr_tt(s-1,:),scale_ratio_k,tgs,s);
				else
					pyr_kernel{s} = downsample_kernel(pyr_kernel{s},pyr_tt(s,:),pyr_kernel{1},pyr_tt(1,:),scale_ratio_k,tgs,s);
				end
			else
				if chain_downsampling
					pyr_kernel{s} = downsample_kernel(pyr_kernel{s}, {pyr_tt{s,1}/scale_ratio_k,pyr_tt{s,2}/scale_ratio_k,pyr_tt{s,3}/scale_ratio_k},pyr_kernel{s-1},pyr_tt(s-1,:),scale_ratio_k,tgs,s);
				else
					pyr_kernel{s} = downsample_kernel(pyr_kernel{s}, {pyr_tt{s,1}/scale_factor,pyr_tt{s,2}/scale_factor,pyr_tt{s,3}/scale_factor},pyr_kernel{1},pyr_tt(1,:),scale_ratio_k,tgs,s);
				end
			end
        end
		if ~sum(pyr_kernel{s}(:)) == 0
			pyr_kernel{s} = pyr_kernel{s}/sum(pyr_kernel{s}(:));
		end
    else
        pyr_kernel{s} = zeros(2*khh_s+1,2*khw_s+1,2*khd_s+1);
    end
    s = s + 1;
end
