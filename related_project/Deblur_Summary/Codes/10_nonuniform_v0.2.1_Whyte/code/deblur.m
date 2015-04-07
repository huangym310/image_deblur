%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'));

if (~exist('load_mask'))
	LOAD_MASK = 0;
else
	LOAD_MASK = 1;
	load_mask = logical(load_mask(:,:,1));
end

clamp_edges_to_zero = 0;
non_uniform = NON_UNIFORM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess observed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure we're using the old version of imresize
% --- new versions align differently, which messes up the 
%     internal camera calibration of the different scales
imver = ver('images');
if str2num(imver.Version) >= 5.4
	imresizefn = @imresize_old;
else
	imresizefn = @imresize;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate blurred image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SYNTHETIC
	% the obs image is really the sharp one if we are synthetic
	% so copy it to true_im
	% true_im_all = obs_im_all;
	% true_grad_all_x = obs_grad_all_x;
	% true_grad_all_y = obs_grad_all_y;
	% true_grad_all = obs_grad_all;
	true_im = obs_im; % the patch
	% true_grad_scale = obs_grad_scale; % gradient pyramid
	% true_grad_scale_x = obs_grad_scale_x; % x-grad pyr
	% true_grad_scale_y = obs_grad_scale_y; % y-grad pyr

	if non_uniform
		% blur true_im
		obs_im = apply_blur_kernel_mex(double(true_im),size(obs_im(:,:,1)),Kblurry_fullres,Kblurry_fullres,-theta_trajectory,blur_kernel(:),clamp_edges_to_zero,non_uniform);
		obs_im = obs_im + SYNTHETIC_NOISE_STD*randn(size(obs_im));
		%%% Now round back to unit8 image
		obs_im = uint8(obs_im);
	else
		if EXTERNAL_SYNTHETIC
			load(fullfile(OUTPUT_DIRECTORY,strrep(OUTPUT_FILENAME,'_uni_','_rot_'),'results.mat'),'obs_im','true_im');
		else
			% blur true_im
			if FFT_SHIFT
				bktmp = ifftshiftpad(blur_kernel,size(true_im));
			else
				bktmp = blur_kernel;
			end
			obs_im = real(ifft2(fft2(true_im) .* fft2(bktmp,size(true_im,1),size(true_im,2))));
			% add noise
			obs_im = obs_im + SYNTHETIC_NOISE_STD*randn(size(obs_im));
			%%% Now round back to unit8 image
			obs_im = uint8(obs_im);
		end
	end
	% Save blurry and original images out
	save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'),'obs_im','true_im','-append','-v7.3');
end

%%% Save original image
obs_im_orig = obs_im;

% Resize image if specified
if PRESCALE && PRESCALE_BEFORE_GRAYSCALE
    obs_im = imresizefn(obs_im,PRESCALE,'bilinear');
	obs_im_orig = imresizefn(obs_im_orig,PRESCALE,'bilinear');
	if LOAD_MASK
		load_mask = imresizefn(load_mask,PRESCALE,'nearest');
	end
end

%%% ROB NOTE: shouldn't rgb2gray be AFTER gamma correction?
% Convert to grayscale
if (size(obs_im,3)>1)
    if SPECIFIC_COLOR_CHANNEL
        obs_im = obs_im(:,:,SPECIFIC_COLOR_CHANNEL);
		if SYNTHETIC
			true_im = true_im(:,:,SPECIFIC_COLOR_CHANNEL);
		end
    else
        obs_im = rgb2gray_rob(obs_im);
		if SYNTHETIC
	        true_im = rgb2gray_rob(true_im);
		end
    end
end

% Resize image if specified
if PRESCALE && ~PRESCALE_BEFORE_GRAYSCALE
    obs_im = imresizefn(obs_im,PRESCALE,'bilinear');
    obs_im_orig = imresizefn(obs_im_orig,PRESCALE,'bilinear');
	if LOAD_MASK
		load_mask = imresizefn(load_mask,PRESCALE,'nearest');
	end
end

% Get size
[obs_imy,obs_imx] = size(obs_im);

% Get maximum distance from principal point
if non_uniform
	max_radius = max(sqrt(([1, obs_imx, 1, obs_imx]-x0).^2 + ([1, 1, obs_imy, obs_imy]-y0).^2));
end

% Gamma correction
if (GAMMA_CORRECTION~=1)
    obs_im = (double(obs_im).^(GAMMA_CORRECTION))/(256^(GAMMA_CORRECTION-1));
	if SYNTHETIC
		true_im = (double(true_im).^(GAMMA_CORRECTION))/(256^(GAMMA_CORRECTION-1));
	end
else
    obs_im = double(obs_im);
	if SYNTHETIC
		true_im = double(true_im);
	end
end

% Rescale image values
if (INTENSITY_SCALING)
    obs_im = obs_im.*INTENSITY_SCALING;
	if SYNTHETIC
		true_im = true_im.*INTENSITY_SCALING;
	end
end

% Find saturated regions of blurry image
if (SATURATION_MASK)
    %%% Use intensity channel
    sat = (obs_im(:,:,1) > 256*(SATURATION_THRESHOLD/256)^GAMMA_CORRECTION);
	if non_uniform
	    q = conv2(double(sat), ones(1+2*[max(blur_y_lims(2),blur_z_lims(2))-min(blur_y_lims(1),blur_z_lims(1)),...
	    							     max(blur_x_lims(2),blur_z_lims(2))-min(blur_x_lims(1),blur_z_lims(1))]),'same');
	else
	    q = conv2(double(sat), ones(size(blur_kernel)),'same');
	end
    mask = (q>0);
else
    mask = zeros(size(obs_im,1),size(obs_im,2));
end

% Optionally load a mask of pixels to use
if (LOAD_MASK)
	mask = mask | ~load_mask;
end

% Images are linear so
if AUTOMATIC_PATCH
    [obs_im_tmp,PATCH_LOCATION] = automatic_patch_selector(obs_im.^(1/DEFAULT_GAMMA),max(PATCH_SIZE),AUTOMATIC_PATCH_CENTER_WEIGHT,mask);
end

if strcmp(CAMERA_TYPE,'unknown')
    %%% leave [0:255] alone
else
    %%% load up map for file
    load(CAMERA_TYPE);
    im_c = obs_im(:);
	new_col = interp1q([0:255]',response_curves(4,:)',im_c);
	obs_im = reshape(new_col,obs_imy,obs_imx);
end


% final level is high res version
obs_im_all = obs_im;
mask_all   = mask;
if SYNTHETIC
	true_im_all = true_im;
end


% select gradient filter
kx = [1 -1]; ky = [1 -1]';
obs_grad_all_x = conv2(obs_im_all,kx,'valid');
obs_grad_all_y = conv2(obs_im_all,ky,'valid');

yy = min(size(obs_grad_all_x,1),size(obs_grad_all_y,1));
xx = min(size(obs_grad_all_x,2),size(obs_grad_all_y,2));

obs_grad_all = [obs_grad_all_x(1:yy,1:xx),obs_grad_all_y(1:yy,1:xx)];
if SYNTHETIC
	true_grad_all_x = conv2(true_im_all,kx,'valid');
	true_grad_all_y = conv2(true_im_all,ky,'valid');
	yy = min(size(true_grad_all_x,1),size(true_grad_all_y,1));
	xx = min(size(true_grad_all_x,2),size(true_grad_all_y,2));
	true_grad_all = [true_grad_all_x(1:yy,1:xx),true_grad_all_y(1:yy,1:xx)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chop out patch and form scale pyramid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Chop out patch in intensity and gradient space
py = 0;
px = 0;
patchrows = PATCH_LOCATION(2)-py:PATCH_LOCATION(2)+PATCH_SIZE(2)-1+py;
patchcols = PATCH_LOCATION(1)-px:PATCH_LOCATION(1)+PATCH_SIZE(1)-1+px;
obs_im = obs_im_all(patchrows,patchcols,:);
if SYNTHETIC
	true_im = true_im_all(patchrows,patchcols,:);
end

mask = mask_all(patchrows,patchcols);

Kblurry = [1, 0, 1-patchcols(1); 0, 1, 1-patchrows(1); 0, 0, 1]*Kblurry_all;

% TODO, should x and y gradients have different intrinsic calibration
% matrices? since e.g. the forward difference method approximates the
% gradient at the midpoint between two adjacent pixels: 
% obs_grad_x(y,x) = obs(y,x+1)-obs(y,x) \approxeq grad_x I (y,x+0.5)
% obs_grad_y(y,x) = obs(y+1,x)-obs(y,x) \approxeq grad_y I (y+0.5,x)

%% chop out gradient patch too
obs_grad_x = obs_grad_all_x(patchrows,patchcols,:);
obs_grad_y = obs_grad_all_y(patchrows,patchcols,:);
if SYNTHETIC
	true_grad_x = true_grad_all_x(patchrows,patchcols,:);
	true_grad_y = true_grad_all_y(patchrows,patchcols,:);
end

if RESCALE_THEN_GRAD
	%%%% rescale images then take gradient

	% use matlab's imresizefn function
	for s = 2:NUM_SCALES
		% resize intrinsic calibration
		Kblurry_scale{NUM_SCALES-s+1} = htranslate([0.5,0.5])*hscale([(1/RESIZE_STEP)^(s-1),(1/RESIZE_STEP)^(s-1)])*htranslate([-0.5,-0.5])*Kblurry;

		obs_im_scale{NUM_SCALES-s+1} = imresizefn(obs_im,(1/RESIZE_STEP)^(s-1),'bilinear');
		if SYNTHETIC
			true_im_scale{NUM_SCALES-s+1} = imresizefn(true_im,(1/RESIZE_STEP)^(s-1),'bilinear');
		end

		%%% Compute saturation mask
		mask_scale{NUM_SCALES-s+1} = ceil(abs(imresizefn(mask,(1/RESIZE_STEP)^(s-1),'nearest')));
		mask_scale{NUM_SCALES-s+1} = mask_scale{NUM_SCALES-s+1}(2:end-1,2:end-1);
	end


	% final level is high res version
	obs_im_scale{NUM_SCALES} = obs_im;
	if SYNTHETIC
		true_im_scale{NUM_SCALES} = true_im;
	end
	mask_scale{NUM_SCALES} = mask(2:end-1,2:end-1);
	Kblurry_scale{NUM_SCALES} = Kblurry;

	% select gradient filter
	kx = [0 1 -1]; ky = [0 1 -1]';

	% apply to resized images
	for s = 1:NUM_SCALES
		obs_grad_scale_x{s} = conv2(obs_im_scale{s},kx,'valid');
		obs_grad_scale_y{s} = conv2(obs_im_scale{s},ky,'valid');
		if SYNTHETIC
			true_grad_scale_x{s} = conv2(true_im_scale{s},kx,'valid');
			true_grad_scale_y{s} = conv2(true_im_scale{s},ky,'valid');
		end
		% cropping changes intrinsic calibration
		Kblurry_grad_scale{s} = [1, 0, -1; 0, 1, -1; 0, 0, 1]*Kblurry_scale{s};
	end

else %%% GRAD THEN RESCALE
	% final level is high res version
	obs_grad_scale_x{NUM_SCALES} = obs_grad_x;
	obs_grad_scale_y{NUM_SCALES} = obs_grad_y;
	if SYNTHETIC
		true_grad_scale_x{NUM_SCALES} = true_grad_x;
		true_grad_scale_y{NUM_SCALES} = true_grad_y;
	end
	mask_scale{NUM_SCALES}  = mask;
	Kblurry_scale{NUM_SCALES} = Kblurry;

	if MANUAL_KERNEL
		manual_kernel{NUM_SCALES} = approximate_kernel; %%% loaded up in script
	end

	% use matlab's imresizefn function
	for s = 2:NUM_SCALES
		% resize intrinsic calibration
		Kblurry_scale{NUM_SCALES-s+1} = htranslate([0.5,0.5])*hscale([(1/RESIZE_STEP)^(s-1),(1/RESIZE_STEP)^(s-1)])*htranslate([-0.5,-0.5])*Kblurry;
		% resize image
		obs_grad_scale_x{NUM_SCALES-s+1} = imresizefn(obs_grad_x,(1/RESIZE_STEP)^(s-1),'bilinear');
		obs_grad_scale_y{NUM_SCALES-s+1} = imresizefn(obs_grad_y,(1/RESIZE_STEP)^(s-1),'bilinear');
		if SYNTHETIC
			true_grad_scale_x{NUM_SCALES-s+1} = imresizefn(true_grad_x,(1/RESIZE_STEP)^(s-1),'bilinear');
			true_grad_scale_y{NUM_SCALES-s+1} = imresizefn(true_grad_y,(1/RESIZE_STEP)^(s-1),'bilinear');
		end
		%%% Compute saturation mask
		mask_scale{NUM_SCALES-s+1} = ceil(abs(imresizefn(mask,(1/RESIZE_STEP)^(s-1),'nearest')));
	end
	% No cropping, so intrinsic calibration of gradient images is same
	% as intensity images
	Kblurry_grad_scale = Kblurry_scale;
end

if RESCALE_THEN_GRAD
	%%% Concatenate the gradient images together
	for s = 1:NUM_SCALES
		obs_grad_scale{s} = [ obs_grad_scale_x{s}(2:end-1,:) , obs_grad_scale_y{s}(:,2:end-1) ];
		if SYNTHETIC
			true_grad_scale{s} = [ true_grad_scale_x{s}(2:end-1,:) , true_grad_scale_y{s}(:,2:end-1) ];
		end
	end
else
	%%% Concatenate the gradient images together
    for s = 1:NUM_SCALES
        obs_grad_scale{s} = [ obs_grad_scale_x{s} , obs_grad_scale_y{s} ];
		if SYNTHETIC
			true_grad_scale{s} = [ true_grad_scale_x{s} , true_grad_scale_y{s} ];
		end
    end
end	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess blur kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if non_uniform
	% Make angle and kernel pyramid
	dm = max(abs([1-x0, obs_imx-x0, 1-y0, obs_imy-y0]));
	% theta_grid_size = calculate_theta_grid_size([obs_imy,obs_imx],Kblurry_all,pixels_per_theta_step);
	theta_grid_size = pixels_per_theta_step * fl / (fl^2 + dm^2);
	theta_stepx = theta_grid_size;
	theta_stepy = theta_grid_size;
	theta_stepz = pixels_per_theta_step / max_radius;
	theta_x_lims = blur_x_lims * fl / (fl^2 + dm^2);
	theta_y_lims = blur_y_lims * fl / (fl^2 + dm^2);
	theta_z_lims = blur_z_lims / max_radius;
	chain_downsampling = true;
	[pyr_kernel,pyr_theta_grid,pyr_theta_step] = make_kernel_pyramid(theta_x_lims,theta_y_lims,theta_z_lims,[theta_stepy,theta_stepx,theta_stepz],1/RESIZE_STEP,NUM_SCALES,1,blur_kernel,chain_downsampling,theta_trajectory);
	% reverse order of pyramid due to different convention, and squash 3D kernel to 1D
	for s=1:length(pyr_kernel)
		blur_kernel_scale{NUM_SCALES-s+1} = pyr_kernel{s}(:);
		theta_scale{NUM_SCALES-s+1} = [pyr_theta_grid{s,2}(:), pyr_theta_grid{s,1}(:), pyr_theta_grid{s,3}(:)]';
		theta_grid_scale(NUM_SCALES-s+1,:) = pyr_theta_grid(s,:);
		theta_step_scale{NUM_SCALES-s+1} = pyr_theta_step{s};
	end
	blur_kernel = blur_kernel_scale{NUM_SCALES};
else
    blur_kernel_scale{NUM_SCALES} = blur_kernel;
    %%% use matlab's imresizefn function
    for s = 2:NUM_SCALES
        dims = size(blur_kernel_scale{NUM_SCALES}) * (1/RESIZE_STEP)^(s-1);
        dims = dims + (1-mod(dims,2)); %% make odd size
        if (min(dims)<4)
            %% Blur blur_kernel first, then resize (since it will use nearest neighbour)
            h = fspecial('gaussian',dims,1);
            blur_kernel_scale{NUM_SCALES-s+1} = imresizefn(conv2(blur_kernel,h),dims,'nearest');
        else
			blur_kernel_scale{NUM_SCALES-s+1} = imresizefn(blur_kernel,dims,'bilinear');
        end
        blur_kernel_scale{NUM_SCALES-s+1} = blur_kernel_scale{NUM_SCALES-s+1} / sum(blur_kernel_scale{NUM_SCALES-s+1}(:));
		% Dummy variables to prevent errors
		theta_step_scale{s} = zeros(1,3);
		theta_grid_scale(s,1:3) = cell(1,3);
    end	
	% Dummy variables to prevent errors
	theta_step_scale{1} = zeros(1,3);
	theta_grid_scale(1,1:3) = cell(1,3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate blurred image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~SYNTHETIC
	%%% NOT SYNTHIETC
	% No idea what this was for... maybe offset image to have correct offsets to use fft for convolution?
	if ~non_uniform
		for s=1:NUM_SCALES
			obs_grad_scale_old{s} = obs_grad_scale{s};
			db = delta_kernel_uni(size(blur_kernel_scale{s},1));
			obs_grad_scale{s} = real(ifft2(fft2(obs_grad_scale{s}) .* ...
										  fft2(db,size(obs_grad_scale{s},1),size(obs_grad_scale{s},2))));
			%%% translate mask too...
			if FERGUS_MASK_TRANSLATION
				% warning('Applying Fergus'' translation of mask for compatibility')
				mask_scale{s} = conv2(double(mask_scale{s}),db,'same');
			else
				% warning('Applying different translation of mask to Fergus''')
				mask_scale{s} = real(ifft2(fft2(mask_scale{s}) .* ...
										   fft2(db,size(mask_scale{s},1),size(mask_scale{s},2))));
			end
		end
	end
end


Ksharp_scale = Kblurry_scale;
Ksharp_grad_scale = Kblurry_grad_scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get sizes of everything and make data and masks for all scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = START_SCALE:NUM_SCALES
    % blur
    K(s) = size(blur_kernel_scale{s},1);
    L(s) = size(blur_kernel_scale{s},2);
	hK(s) = floor(K(s)/2); hL(s) = floor(L(s)/2);

    % image
    M(s) = size(obs_grad_scale{s},1);
    N(s) = size(obs_grad_scale{s},2);

    if BLUR_MASK
		if non_uniform
			error('Oliver didnt implement this')
		end
		for a=1:BLUR_COMPONENTS
			[xx,yy] = meshgrid(-hK(s):hK(s),-hL(s):hL(s));
			spatial_blur_mask{s}(a,:) = log(normMDpdf([xx(:)';yy(:)'],[0 0]',eye(2)*K(s)*BLUR_MASK_VARIANCES(a)));
		end
    else
        spatial_blur_mask{s} = zeros(BLUR_COMPONENTS,K(s)*L(s));
    end

	spatial_image_mask{s} = [zeros(size(mask_scale{s})),zeros(size(mask_scale{s}))];

    % observed
    if ~non_uniform
		I(s) = 2*M(s);
        J(s) = 2*N(s);
        Dp{s} = zeros(I(s),J(s));
		if (FFT_SHIFT)
			Dp{s}(hK(s)+1:M(s)-hK(s),hL(s)+1:N(s)/2-hL(s),:) = 1;
			Dp{s}(hK(s)+1:M(s)-hK(s),N(s)/2+hL(s)+1:N(s)-hL(s),:) = 1;
		else
			Dp{s}(K(s):M(s),L(s):N(s)/2) = 1; %% x-plane
			Dp{s}(K(s):M(s),L(s)+N(s)/2:N(s)) = 1; %% y-plane
		end

        D{s} = padarray(obs_grad_scale{s},[M(s) N(s)],0,'post');

        %%% Now add in saturation mask
        if (SATURATION_MASK==1)
            Dp{s} = Dp{s} .* padarray(1-[mask_scale{s},mask_scale{s}],[M(s) N(s)],0,'post');
        end
    else
        I(s) = M(s);
        J(s) = N(s);

        D{s} = obs_grad_scale{s};
        % make masks
		Dp{s} = zeros(I(s),J(s));
		Dp0{s} = mask_observation_border([I(s),J(s)/2],Kblurry_grad_scale{s},theta_grid_scale{s,2},theta_grid_scale{s,1},theta_grid_scale{s,3},Ksharp_scale{s},[I(s),J(s)/2]);
		% Dp0{s} = apply_blur_kernel_mex(ones([I(s),J(s)/2]),[I(s),J(s)/2],Kblurry_grad_scale{s},Kblurry_grad_scale{s},-theta_scale{s},ones(1,K(s)),clamp_edges_to_zero,non_uniform);
		% Dp0{s} = mask_observation_border_mex(tmp,[I(s),J(s)/2],Kblurry_grad_scale{s},Kblurry_grad_scale{s},-theta_scale{s},ones(1,K(s)),1,non_uniform);
		if max(Dp0{s}(:)) > 1 % if using mask_observation_border_mex
			Dp{s}(1:M(s),1:N(s)/2) = Dp0{s} > K(s)*L(s)-0.5;
			Dp{s}(1:M(s),N(s)/2+1:N(s)) = Dp0{s} > K(s)*L(s)-0.5;
		else
			Dp{s}(1:M(s),1:N(s)/2) = Dp0{s};
			Dp{s}(1:M(s),N(s)/2+1:N(s)) = Dp0{s};
		end
		%%% Now add in saturation mask
		if (SATURATION_MASK==1)
			Dp{s} = Dp{s} .* (1-[mask_scale{s},mask_scale{s}]);
		end
    end
end

if RESUME_PARTIAL
	load(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results_partial.mat'));
	if exist('START_SCALE_RESUME')
		START_SCALE = START_SCALE_RESUME;
	else
		START_SCALE = length(me_est) + 1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_ESTIMATE_KERNEL
	for s = START_SCALE:NUM_SCALES
	    if ((NUM_SCALES-s+1)>8)
	        priors(NUM_SCALES-s+1) = priors(8);
	    end

	    % Loop over scales
	    %             nDims,size_per_dims,#prior compoents,#prior type,lock prior,update]
	    dimensions = [1       1         1                0           0         1;
	        		  1       K(s)*L(s) BLUR_COMPONENTS  BLUR_PRIOR  BLUR_LOCK 1;
	        		  1		  M(s)*N(s) IMAGE_COMPONENTS IMAGE_PRIOR 1         1];

		options = struct('converge_criteria',CONVERGENCE, ...
						 'plot_linesearch',PLOT_LINESEARCH, ...
						 'DISPLAY_EACH_ITERATION',DISPLAY_EACH_ITERATION, ...
						 'SAVE_EACH_ITERATION',SAVE_EACH_ITERATION, ...
						 'init_noise_var',NOISE_INIT, ...
						 'restart_priors',0, ...
						 'restart_switched_off',0, ...
						 'num_iters',MAX_ITERATIONS, ...
						 'text',['Scale=' int2str(s) '/' int2str(NUM_SCALES)], ...
						 'non_uniform',non_uniform, ...
	   					 'I',I(s), ...
	   					 'J',J(s), ...
	   					 'K',K(s), ...
	   					 'L',L(s), ...
	   					 'M',M(s), ...
	   					 'N',N(s), ...
	   					 'priors',priors(NUM_SCALES-s+1), ...
						 'FFT_SHIFT',FFT_SHIFT, ...
						 'nLayers',1, ...
						 'theta_step',theta_step_scale{s}, ...
						 'Kblurry',Kblurry_grad_scale{s}, ...
						 'Ksharp',Kblurry_grad_scale{s}, ...
						 'OUTPUT_DIRECTORY',OUTPUT_DIRECTORY, ...
						 'OUTPUT_FILENAME',OUTPUT_FILENAME, ...
						 's',s, ...
						 'theta_grid',{theta_grid_scale(s,:)}, ...
						 'NUM_THREADS',NUM_THREADS);
	    % Initialize
	    if (s==START_SCALE && ~RESUME_PARTIAL) % first iteration
			initoptions = struct('pres',INIT_PRESCISION,...
								 'prior',IMAGE_PRIOR,...
								 'prior_num',IMAGE_COMPONENTS,...
								 'mode_im',FIRST_INIT_MODE_IMAGE,...
								 'mode_blur',FIRST_INIT_MODE_BLUR);
								
	        if SYNTHETIC
				if non_uniform
					[x1,x2] = initialize_parameters2_rot(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s}, true_grad_scale{s},initoptions,spatial_image_mask{s},options);
				else
		            [x1,x2] = initialize_parameters2_uni(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s}, true_grad_scale{s},initoptions,spatial_image_mask{s},options);
				end
	        else
				if non_uniform
					[x1,x2] = initialize_parameters2_rot(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s}, obs_grad_scale{s},initoptions,spatial_image_mask{s},options);
				else
		            [x1,x2] = initialize_parameters2_uni(obs_grad_scale{s},blur_kernel_scale{s},obs_grad_scale{s},blur_kernel_scale{s}, obs_grad_scale{s},initoptions,spatial_image_mask{s},options);
				end
	        end
			use_rotations_scale{s} = true(K(s)*L(s),1);
	    else % as we work up pyramid
			initoptions = struct('pres',INIT_PRESCISION,...
								 'prior',IMAGE_PRIOR,...
								 'prior_num',IMAGE_COMPONENTS,...
								 'mode_im',INIT_MODE_IMAGE,...
								 'mode_blur',INIT_MODE_BLUR);
	        if SYNTHETIC
				if non_uniform
					[x1,x2] = initialize_parameters2_rot(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},blur_kernel_scale{s}, true_grad_scale{s},initoptions,spatial_image_mask{s},options);
				else
		            [x1,x2] = initialize_parameters2_uni(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},blur_kernel_scale{s}, true_grad_scale{s},initoptions,spatial_image_mask{s},options);
				end
	        else
				if non_uniform
					[x1,x2] = initialize_parameters2_rot(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},[],[],initoptions,spatial_image_mask{s},options);
				else
		            [x1,x2] = initialize_parameters2_uni(obs_grad_scale{s},new_blur{s-1},new_grad{s-1},[],[],initoptions,spatial_image_mask{s},options);
				end
	        end
	    end

	    m = x1./x2;
	    mx_init{s} = reshape(m(2+K(s)*L(s):end),M(s),N(s));
	    me_init{s} = reshape(m(2:1+K(s)*L(s)),K(s),L(s));

	    % Call routine
	    [ensemble,D_log{s},gamma_log] = train_ensemble_main6(dimensions,x1,x2,options,D{s},Dp{s},spatial_blur_mask{s},1-(spatial_image_mask{s}(:)>0),use_rotations_scale{s});

	    % Extract mean blur/image
	    me_est{s} =reshape(train_ensemble_get(2,dimensions,ensemble.mx),K(s),L(s));
	    e1_est{s} =reshape(train_ensemble_get(2,dimensions,ensemble.x1),K(s),L(s));
	    e2_est{s} =reshape(train_ensemble_get(2,dimensions,ensemble.x2),K(s),L(s));

	    mx_est{s} =reshape(train_ensemble_get(3,dimensions,ensemble.mx),M(s),N(s));
	    x1_est{s} =reshape(train_ensemble_get(3,dimensions,ensemble.x1),M(s),N(s));
	    x2_est{s} =reshape(train_ensemble_get(3,dimensions,ensemble.x2),M(s),N(s));

		ensemble_est{s} = ensemble;

	    if (s~=NUM_SCALES)
			me_est_tmp = me_est{s};
			if (THRESHOLD_SCALE_KERNEL)
				me_est_tmp(me_est_tmp < max(me_est_tmp(:))/KERNEL_THRESHOLD) = 0;
			end

	        % Use solution for next level up...
			moveoptions = struct('K',K(s+1), ...
								 'L',L(s+1), ...
								 'M',M(s+1), ...
								 'N',N(s+1), ...
								 'RESIZE_STEP',RESIZE_STEP, ...
								 'CENTER_BLUR',CENTER_BLUR);
			if non_uniform
				[new_grad{s},new_blur{s}] = move_level_rot(mx_est{s},me_est_tmp,moveoptions,theta_step_scale{s},Ksharp_scale{s},theta_grid_scale(s,:),theta_grid_scale(s+1,:));
				new_grad{s} = new_grad{s}*numel(e1_est{s})/numel(new_grad{s});
			else
				[new_grad{s},new_blur{s}] = move_level_uni(mx_est{s},me_est_tmp,moveoptions);
			end

	        if MANUAL_KERNEL
	            % Overwrite upsampled version with manual version
	            new_blur{s} = blur_kernel_scale{s+1};
	        end

			if MASK_KERNEL
				if non_uniform
					use_rotations_scale{s+1} = reshape(new_blur{s},size(theta_grid_scale{s+1,1})) > 0;
				else
					use_rotations_scale{s+1} = reshape(new_blur{s},[K(s+1),L(s+1)]) > 0;
				end
	            nkdims = length(size(use_rotations_scale{s+1}));
	            if KERNEL_DILATE_RADIUS > 0
	                if KERNEL_DILATE_RADIUS < inf
	                    dilate_kernel = ones([repmat(1+2*KERNEL_DILATE_RADIUS,[1, nkdims]), 1]);
	                    use_rotations_scale{s+1} = convn(double(use_rotations_scale{s+1}),dilate_kernel,'same') >= 0.9;
	                else
						if non_uniform
							use_rotations_scale{s+1} = true(numel(theta_grid_scale{s+1,1}),1);
						else
							use_rotations_scale{s+1} = true(K(s+1),L(s+1));
	                    end
	                end
	            end
			else
				if non_uniform
					use_rotations_scale{s+1} = true(numel(theta_grid_scale{s+1,1}),1);
				else
					use_rotations_scale{s+1} = true(K(s+1),L(s+1));
	            end
			end
			if non_uniform
				use_rotations_scale{s+1} = logical(use_rotations_scale{s+1}(:));
			else
				use_rotations_scale{s+1} = logical(use_rotations_scale{s+1});
			end			

	    end

	    %%% save everything out
		save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results_tmp.mat'), 'mx_init', 'me_init', 'mx_est', 'me_est', 'new_grad', 'new_blur', 'ensemble', 'D_log', 'theta_grid_scale', 'D', 'Dp', 'Kblurry_scale', 'use_rotations_scale', 'dimensions', 'ensemble_est', '-v7.3');

		figure(1000+s);
		subplot(121); imagesc(reconsEdge3(mx_est{s})); axis image; axis off; colormap gray; title(['Image, scale ',num2str(s)]);
		if non_uniform
			subplot(122); plot_nonuni_kernel(me_est{s},theta_grid_scale(s,:),1,0,0,0,1); title(['Kernel, scale ',num2str(s)]);
		else
			subplot(122); imagesc(me_est{s}); axis image; axis off; colormap gray; title(['Kernel, scale ',num2str(s)]);
		end
		drawnow
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Reconstruct image patch to intensity space from gradients
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	obs_im_recon = reconsEdge3(mx_est{NUM_SCALES});

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save to mat file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'),...
	    'ensemble',...
	    'dimensions',...
	    'obs_im_recon',...
	    'obs_im_orig',...
	    'obs_im_all',...
	    'obs_grad_scale',...
	    'mask_all',...
	    'blur_kernel_scale',...
	    'use_rotations_scale',...
	    'obs_im_recon',...
	    'blur_kernel',...
	    'mx_init',...
	    'me_init',...
	    'mx_est',...
	    'me_est',...
	    'new_grad',...
	    'new_blur',...
	    'PATCH_LOCATION',...
	    'theta_grid_scale',...
	    '-append','-v7.3');
else
	load(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results_tmp.mat'));
	load(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'));
	if ~exist('use_rotations_scale')
		for ss=1:length(me_est)
			use_rotations_scale{ss} = true(size(me_est{ss}));
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now do actual deblurring of whole image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_DECONVOLUTION
	if strcmp(IMAGE_RECONSTRUCTION,'lucy')
	    %%% Run RL script
		if non_uniform
			[deblurred_im_lucy,kernel_out_lucy] = fiddle_lucy3_rot(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME),LUCY_ITS,1,SCALE_OFFSET,KERNEL_THRESHOLD);
		else
			[deblurred_im_lucy,kernel_out_lucy] = fiddle_lucy3_uni(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME),LUCY_ITS,1,SCALE_OFFSET,KERNEL_THRESHOLD);
		end
	else
	    error('Unrecognised type of deblurring');
	end

	% save final image and kernel
	if strcmp(IMAGE_RECONSTRUCTION,'lucy')
		save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'),'deblurred_im_lucy','kernel_out_lucy','-append');
	end
end
