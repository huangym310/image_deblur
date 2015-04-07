function i_rl = deconvlucy_rotational(imblur,dimsout,blur_kernel,theta_list,Kblurry,Ksharp,sat_thresh,num_iters)
% 	deconvlucy_rotational		Apply Richardson-Lucy deconvolution within the rotational camera shake model
% 		i_rl = deconvlucy_rotational(imblur,dimsout,blur_kernel,theta_list,Kblurry,Ksharp,sat_thresh,num_iters)
% 
% 		Inputs:
% 				imblur			Blurry image
% 				dimsout			[height,width] of sharp image
% 				blur_kernel		Blur kernel to deconvolve with
% 				theta_list		3 x numel(blur_kernel) array, containing orientations to which elements of blur_kernel correspond
% 				Kblurry			Internal calibration matrix of blurry image
% 				Ksharp			Internal calibration matrix of sharp image
% 				sat_thresh		Threshold above which, a pixel will be considered saturated
% 				num_iters		Number of Richarson-Lucy iterations
% 
% 		Outputs:
% 				i_rl			Deblurred image
% 
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

if nargin < 8, num_iters = 10; end
heio = dimsout(1);
wido = dimsout(2);
[heis,wids,channels] = size(imblur);

% Remove the zero elements of kernel
use_rotations = blur_kernel(:) ~= 0;
blur_kernel = blur_kernel(use_rotations);
theta_list = theta_list(:,use_rotations);

% Calculate padding
imcornersshake = [1,   1,wids,wids;...
                  1,heis,heis,   1;...
                  1,   1,   1,   1];
hfs_x1 = 0;
hfs_x2 = 0;
hfs_y1 = 0;
hfs_y2 = 0;
% for each non-zero in the kernel...
for i=1:size(theta_list,2)
    % back proect corners of shaken image to see how far out we need to pad
    H = Ksharp*expm(crossmatrix(-theta_list(:,i)))*inv(Kblurry);
    imcornersnoisy = hnormalise(H*imcornersshake);
    offsets = abs(imcornersnoisy-imcornersshake);
    if offsets(1,1) > hfs_x1, hfs_x1 = ceil(offsets(1,1)); end
    if offsets(1,2) > hfs_x1, hfs_x1 = ceil(offsets(1,2)); end
    if offsets(1,3) > hfs_x2, hfs_x2 = ceil(offsets(1,3)); end
    if offsets(1,4) > hfs_x2, hfs_x2 = ceil(offsets(1,4)); end
    if offsets(2,1) > hfs_y1, hfs_y1 = ceil(offsets(2,1)); end
    if offsets(2,2) > hfs_y2, hfs_y2 = ceil(offsets(2,2)); end
    if offsets(2,3) > hfs_y1, hfs_y1 = ceil(offsets(2,3)); end
    if offsets(2,4) > hfs_y2, hfs_y2 = ceil(offsets(2,4)); end
end

Kblurry = [1, 0, hfs_x1;...
     0, 1, hfs_y1;...
     0, 0, 1     ]...
     * Kblurry;
 
Ksharp = [1, 0, hfs_x1;...
     0, 1, hfs_y1;...
     0, 0, 1     ]...
     * Ksharp;

wids=wids+hfs_x1+hfs_x2;
heis=heis+hfs_y1+hfs_y2;
wido=wido+hfs_x1+hfs_x2;
heio=heio+hfs_y1+hfs_y2;

imblur = imblur([ones(1,hfs_y1),1:end,end*ones(1,hfs_y2)],[ones(1,hfs_x1),1:end,end*ones(1,hfs_x2)],:);

theta_pre = [0;0;0];
clamp_edges_to_zero = 0;
non_uniform = 1;

blurfn     = @(im) apply_blur_kernel_mex(im,[heis,wids],Ksharp,Kblurry,-theta_list,blur_kernel,clamp_edges_to_zero,non_uniform,-theta_pre);
blurconjfn = @(im) apply_blur_kernel_mex(im,[heio,wido],Kblurry,Ksharp, theta_list,blur_kernel,clamp_edges_to_zero,non_uniform,-theta_pre);

[tmp,center_csf_element] = min(sum(theta_list.^2,1));
tmpcsf_reproject = zeros(size(blur_kernel));
tmpcsf_reproject(center_csf_element) = 1;
% Initialise sharp image by projecting blurry image into sharp frame
i_rl = apply_blur_kernel_mex(imblur,[heio,wido],Kblurry,Ksharp, theta_list,tmpcsf_reproject,clamp_edges_to_zero,non_uniform,-theta_pre);

for iter = 1:num_iters
    i_rl = blurconjfn(imblur./max(blurfn(i_rl),eps)) .* i_rl;
end

i_rl = i_rl(hfs_y1+1:end-hfs_y2,hfs_x1+1:end-hfs_x2,:);

