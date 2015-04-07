function [out,blur_kernel]=fiddle_lucy3_uni(file_name,LUCY_ITS_IN,SAVE_TO_DISK,SCALE_OFFSET_IN,THRESHOLD_IN)
%
% Rountine to call Richardson-Lucy algorithm after kernel inference  
%
% Inputs:
% LUCY_ITS - integer. Default = 10. Sensible values 5 to
%    30. Number of Lucy-Richardson iterations to use when unblurring
%    the image using the inferred kernel. 10 is the default but if the
%    kernel is very long and thin you might need more. If you turn it
%    up too much it amplifies the noise in the image.
%
% SAVE_TO_DISK - binary. Turns saving of deblurred image and final kernel
% (post thresholding) to disk or not.
%
%
% SCALE_OFFSET - integer. Default = 0. Sometimes it may not be
%      possible to deblur the image at full resolution due to high
%      noise levels (see discussion in 7.3). If this is the case, you
%      can tell the fiddle_lucy3 function to use a coarser scale. The
%      value of SCALE_OFFSET dicates how many scale levels you drop
%      down before selecting the kernel. i.e. SCALE_OFFSET = 1 will use
%      a kernel that is one scale level (RESIZE_STEP smaller) than the
%      full resolution kernel.
%
% THRESHOLD - float. Default = 7, sensible range is 5 to
%    15. This is a threshold on kernel intensities, applied after the
%    whole inference stage, designed to remove noise from the
%    kernel. It is a dynamic thresold which is the percentage of the
%    maximum value of the kernel. Since it is a bit dependent on the
%    intensity profile of the kernel, some manual tuning may be needed
%    to get the optimal value. If you don't want to use the threshold at all, set it to a
%    large number (i.e. 100). If you set it too high, it will start to
%    erode the structure of the kernel. This parameter is a bit of a
%    hack - in reality you shouldn't need it, but it can make quite a
%    bit of difference if the inferred kernel is noisy. 
%
% Outputs:
% 1. out - Deblurred image  
% 2. blur_kernel - blur kernel after thresholding step.

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% Modified by Oliver Whyte for CVPR 2010 paper:
% "Non-uniform Deblurring for Shaken Images"
% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

imver = ver('images');
if str2num(imver.Version) >= 5.4
	imresizefn = @imresize_old;
else
	imresizefn = @imresize;
end
  
SHOW_BOX = 0; % show gray rectangle on blurry image indicating selected region
EDGE_CROP = 0; % crop edge of deblurred image off (since you can't
               % recover it).

			% obs_jpg = imread('~/Work/deblur/fergus/images/newsofa_raw.jpg'); warning('loading sofa image');

if (nargin==4)
  THRESHOLD_IN = 10;
end

if (nargin==3)
  SCALE_OFFSET_IN = 0;
  THRESHOLD_IN = 10;
end

if (nargin==2)
  SCALE_OFFSET_IN = 0;
  SAVE_TO_DISK_IN = 0;
  THRESHOLD_IN = 10;
end

if (nargin==1)
  LUCY_ITS_IN = 10;
  SCALE_OFFSET_IN = 0;
  SAVE_TO_DISK_IN = 0;
  THRESHOLD_IN = 10;
end  

%%% Load up model
load(fullfile(file_name,'results.mat'));
if ~exist('mx_est')
	[pp,nn,xx,vv] = fileparts(file_name);
	% load(fullfile(pp,['tmp_',nn,xx,vv]));
	load(fullfile(file_name,'results_tmp.mat'));
	% add empty cells to end of mx_est and me_est, and set SCALE_OFFSET
	ll = length(theta_grid_scale);
	sl = length(mx_est);
	mx_est(sl+1:ll) = {0};
	me_est(sl+1:ll) = {0};
	SCALE_OFFSET_IN = ll - sl;
end

%% override values in save file
SCALE_OFFSET = SCALE_OFFSET_IN;
LUCY_ITS = LUCY_ITS_IN;
THRESHOLD = THRESHOLD_IN;

s = length(theta_grid_scale) - SCALE_OFFSET;

%%% Get blur kernel
blur_kernel = me_est{end-SCALE_OFFSET}/sum(me_est{end-SCALE_OFFSET}(:));

%%% Threshold kernel
threshold = max(blur_kernel(:))/THRESHOLD;
z=find(blur_kernel(:)<threshold);
blur_kernel(z)=0;

% check kernel sums to 1
blur_kernel = blur_kernel / sum(blur_kernel(:));

obs_im = double(obs_im);

%%% Resize image
if (PRESCALE~=1)
  obs_im = imresizefn(obs_im,PRESCALE,'bilinear');
end

%%% Rescale if not using final scale
if (SCALE_OFFSET>0)
  obs_im = imresizefn(obs_im,(1/sqrt(2))^SCALE_OFFSET,'bilinear');
end

if max(obs_im(:)) <= 1
	obs_im = obs_im * 256;
end

%%% Take out gamma corection
if (GAMMA_CORRECTION~=1)
	obs_im_gam = ((double(obs_im).^GAMMA_CORRECTION)/(256^(GAMMA_CORRECTION-1)));
else
	obs_im_gam = double(obs_im);
end
SATURATION_THRESHOLD = (SATURATION_THRESHOLD^GAMMA_CORRECTION)/(256^(GAMMA_CORRECTION-1));

fprintf('.');

% use edgetaper to sort out edge effects
% obs_im_gam = edgetaper(obs_im_gam,blur_kernel);
obs_im_taper = obs_im_gam;
pad_taper = [0,0];
obs_mask = zeros(size(obs_im_gam(:,:,1)));
obs_mask(ceil((size(blur_kernel,1)-1)/2)+1:end-floor((size(blur_kernel,1)-1)/2),ceil((size(blur_kernel,2)-1)/2)+1:end-floor((size(blur_kernel,2)-1)/2)) = 1;

[obs_im_taper,pad_taper] = edgetaper_uni(obs_im_gam,obs_mask,blur_kernel);

fprintf('.');

% run RL
% out = deconvlucy(obs_im_gam,blur_kernel,LUCY_ITS);  
out = deconvlucy_uni(obs_im_taper,blur_kernel,SATURATION_THRESHOLD,LUCY_ITS);


fprintf('.');

dotruedeblur = false;
if ~isempty(strfind(file_name,'levin'))
	dotruedeblur = true;
	levingt = load(fullfile('/path/to/LevinEtalCVPR09Data',[file_name(strfind(file_name,'levin')+6:strfind(file_name,'_uni')-1)]));
	blur_kernel_true = rot90(rot90(levingt.f));
	blur_kernel_true = blur_kernel_true/sum(blur_kernel_true(:));
end

if dotruedeblur
	outtrue = deconvlucy_uni(obs_im_taper,blur_kernel_true,SATURATION_THRESHOLD,LUCY_ITS);
end

     
if (GAMMA_CORRECTION~=1)
  %%% Put gamma back in
  out = 256*((double(out)/256).^(1/GAMMA_CORRECTION));
  if dotruedeblur, outtrue = 256*((double(outtrue)/256).^(1/GAMMA_CORRECTION)); end
else
  out = double(out);
  if dotruedeblur, outtrue = double(outtrue); end
end

% crop from edgetaper
out = out(pad_taper(1)+1:end-pad_taper(1),pad_taper(2)+1:end-pad_taper(2),:);
if dotruedeblur
	outtrue = outtrue(pad_taper(1)+1:end-pad_taper(1),pad_taper(2)+1:end-pad_taper(2),:);
end

for c=1:size(out,3)
	%% shift and rescale to make 0 to 1 image
	out(:,:,c) = out(:,:,c) - min(min(out(:,:,c)));
	out(:,:,c) = out(:,:,c) / max(max(out(:,:,c)));

	%% now do histogram equalization
	if exist('obs_jpg')
		out(:,:,c) = histmatch(out(:,:,c),uint8(imresizefn(obs_jpg(:,:,c),size(out(:,:,1)),'bilinear')));
	else
		out(:,:,c) = histmatch(out(:,:,c),uint8(obs_im(:,:,c)));
	end

	if dotruedeblur
		outtrue(:,:,c) = outtrue(:,:,c) - min(min(outtrue(:,:,c)));
		outtrue(:,:,c) = outtrue(:,:,c) / max(max(outtrue(:,:,c)));
		if exist('obs_jpg')
			outtrue(:,:,c) = histmatch(outtrue(:,:,c),uint8(imresizefn(obs_jpg(:,:,c),size(out(:,:,1)),'bilinear')));
		else
			outtrue(:,:,c) = histmatch(outtrue(:,:,c),uint8(obs_im(:,:,c)));
		end
	end
end

out = uint8(out);
if dotruedeblur
	outtrue = uint8(outtrue);
end

%% apply box to im_col
if SHOW_BOX
  aa = round(AXIS * (1/sqrt(2))^SCALE_OFFSET);
  obs_im(aa(3),aa(1):aa(2),:) = 75;
  obs_im(aa(4),aa(1):aa(2),:) = 75;
  obs_im(aa(3):aa(4),aa(1),:) = 75;
  obs_im(aa(3):aa(4),aa(2),:) = 75;
end

%%% Crop image to avoid edge artifacts
if EDGE_CROP
  edge_offset = floor(size(blur_kernel,1)/2);
  out = out(edge_offset+1:end-edge_offset-1,edge_offset+1:end-edge_offset-1,:);
  obs_im = obs_im(edge_offset+1:end-edge_offset-1,edge_offset+1:end-edge_offset-1,:);
end

%%%%% Show images/kernel

if exist('obs_jpg')
	obs_show = imresizefn(obs_jpg,size(obs_im(:,:,1)),'bilinear');
else
	obs_show = obs_im;
end

%% Plot output image
figure; imagesc(out); title('Output'); axis equal; axis off; colormap gray
%% Plot original, blurry image
figure; imagesc(uint8(obs_show)); title('Original image'); axis equal; axis off; colormap gray
if dotruedeblur
	figure; imagesc(outtrue); title('Deblurred with true kernel'); axis equal; axis off; colormap gray
end

%%% Plot kernels
h=figure; imagesc(blur_kernel); colormap gray; axis square; axis off;

if SAVE_TO_DISK
  % save blurry
  imwrite(uint8(obs_show),fullfile(file_name,sprintf('blurry_rl_s%02d.jpg',s)),'jpg','Quality',100);
  % save blurry
  imwrite(uint8(out),fullfile(file_name,sprintf('deblurred_rl_s%02d.jpg',s)),'jpg','Quality',100);
	if dotruedeblur
		imwrite(uint8(outtrue),fullfile(file_name,sprintf('deblurredtruekernel_rl_s%02d.jpg',s)),'jpg','Quality',100);
	end
  % save kernel
	print(h,'-dpng','-r100',fullfile(file_name,sprintf('kernel_rl_s%02d.png',s)));
	if dotruedeblur
		figure(h); imagesc(blur_kernel_true); axis square; axis off
		print(h,'-dpng','-r100',fullfile(file_name,sprintf('truekernel_rl_s%02d.png',s)));
	end
end
