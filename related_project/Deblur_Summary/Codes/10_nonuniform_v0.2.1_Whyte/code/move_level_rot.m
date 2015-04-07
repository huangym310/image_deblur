function [mx_new,me_new] = move_level_rot(mx,me,moveoptions,theta_step_s,Ksharp,tt_s,tt_up)
	% Author: Rob Fergus
	% Version: 1.0, distribution code.
	% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
	% Copyright 2006, Massachusetts Institute of Technology

	% Modified by Oliver Whyte for CVPR 2010 paper:
	% "Non-uniform Deblurring for Shaken Images"
	% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

	K           = moveoptions.K;
	L           = moveoptions.L;
	M           = moveoptions.M;
	N           = moveoptions.N;
	resize_step = moveoptions.RESIZE_STEP;
	center      = moveoptions.CENTER_BLUR;

	imver = ver('images');
	if str2num(imver.Version) >= 5.4
		imresizefn = @imresize_old;
	else
		imresizefn = @imresize;
	end

	%%% Get number of image channels
	mxc = size(mx,3);

	if (center)
		fprintf('Centering kernel\n');
		me = me / sum(me(:));

		tsy_s = numel(unique(tt_s{1}(:)));
		tsx_s = numel(unique(tt_s{2}(:)));
		tsz_s = numel(unique(tt_s{3}(:)));
		k_s = reshape(me,[tsy_s,tsx_s,tsz_s]);

		k_s0 = k_s;
		mx0 = mx;

		%% get centre of mass
		mu_y = sum((1:tsy_s)' .* reshape(sum(sum(k_s,2),3),[],1));
		mu_x = sum((1:tsx_s)' .* reshape(sum(sum(k_s,1),3),[],1));    
		mu_z = sum((1:tsz_s)' .* reshape(sum(sum(k_s,1),2),[],1));    

		%% get mean offset
		offset_x = round( floor(tsx_s/2)+1 - mu_x );
		offset_y = round( floor(tsy_s/2)+1 - mu_y );
		offset_z = round( floor(tsz_s/2)+1 - mu_z );

		%% make kernel to do translation
		%     shift_kernel = zeros(abs(offset_y*2)+1,abs(offset_x*2)+1,abs(offset_z*2)+1);
		shift_kernel = zeros(tsy_s,tsx_s,tsz_s);
		shift_kernel((tsy_s-1)/2+1+offset_y,(tsx_s-1)/2+1+offset_x,(tsz_s-1)/2+1+offset_z) = 1;

		%% shift both image and blur kernel
		%     me = conv2(me,shift_kernel,'same');
		k_s = convn(k_s,shift_kernel,'same');
		me = reshape(k_s,size(me));

		theta_list = [tt_s{2}(:),tt_s{1}(:),tt_s{3}(:)]';
		clamp_edges_to_zero = 0;
		non_uniform = 1;


		for c=1:mxc
			%       mx(:,:,c) = conv2(mx(:,:,c),flipud(fliplr(shift_kernel)),'same');
			mx(:,1:end/2,c) = apply_blur_kernel_mex(mx(:,1:end/2,c),size(mx(:,1:end/2,c)),Ksharp,Ksharp,theta_list,shift_kernel,clamp_edges_to_zero,non_uniform);
			mx(:,end/2+(1:end/2),c) = apply_blur_kernel_mex(mx(:,end/2+(1:end/2),c),size(mx(:,end/2+(1:end/2),c)),Ksharp,Ksharp,theta_list,shift_kernel,clamp_edges_to_zero,non_uniform);
		end
	end

	mx_new = imresizefn(mx,[M N],'bilinear');
	%     me_new = imresizefn(me,[K L],'bilinear');
	tsy_s = numel(unique(tt_s{1}(:)));
	tsx_s = numel(unique(tt_s{2}(:)));
	tsz_s = numel(unique(tt_s{3}(:)));
	k_s = reshape(me,[tsy_s,tsx_s,tsz_s]);
	tsy_up = numel(unique(tt_up{1}(:)));
	tsx_up = numel(unique(tt_up{2}(:)));
	tsz_up = numel(unique(tt_up{3}(:)));
	k_up = zeros(tsy_up,tsx_up,tsz_up);
	k_up = upsample_kernel(k_s,tt_s,k_up,tt_up,theta_step_s);
	me_new = reshape(k_up,[K L]);

	%%% ensure blur kernel is normalized
	mfactor = sum(me_new(:));
	me_new = me_new / mfactor;

	%%% now 
	%mx_new = mx_new * mfactor;
