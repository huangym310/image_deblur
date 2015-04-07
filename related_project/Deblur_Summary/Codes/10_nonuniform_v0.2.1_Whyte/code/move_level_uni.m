function [mx_new,me_new]=move_level_uni(mx,me,moveoptions)
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
		%% get centre of mass
		mu_y = sum([1:size(me,1)] .* sum(me,2)');
		mu_x = sum([1:size(me,2)] .* sum(me,1));    

		%% get mean offset
		offset_x = round( floor(size(me,2)/2)+1 - mu_x );
		offset_y = round( floor(size(me,1)/2)+1 - mu_y );

		%% make kernel to do translation
		shift_kernel = zeros(abs(offset_y*2)+1,abs(offset_x*2)+1);
		shift_kernel(abs(offset_y)+1+offset_y,abs(offset_x)+1+offset_x) = 1;

		%% shift both image and blur kernel
		me = conv2(me,shift_kernel,'same');

		for c=1:mxc
			mx(:,:,c) = conv2(mx(:,:,c),flipud(fliplr(shift_kernel)),'same');
		end

	end

	mx_new = imresizefn(mx,[M N],'bilinear');
	me_new = imresizefn(me,[K L],'bilinear');   

	%%% ensure blur kernel is normalized
	mfactor = sum(me_new(:));
	me_new = me_new / mfactor;

	%%% now 
	%mx_new = mx_new * mfactor;
