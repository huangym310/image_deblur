function [x1,x2] = initialize_parameters2_uni(obs,blur,im,true_blur,true_im,initoptions,spatial_mask,options)
	% Author: Rob Fergus
	% Version: 1.0, distribution code.
	% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
	% Copyright 2006, Massachusetts Institute of Technology

	% Modified by Oliver Whyte for CVPR 2010 paper:
	% "Non-uniform Deblurring for Shaken Images"
	% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce
	
	pres        = initoptions.pres;
	prior       = initoptions.prior;
	prior_num   = initoptions.prior_num;
	mode_im     = initoptions.mode_im;
	mode_blur   = initoptions.mode_blur;

	priors      = options.priors;
	nLayers     = options.nLayers;
	FFT_SHIFT   = options.FFT_SHIFT;
	Kblurry     = options.Kblurry;
	theta_grid  = options.theta_grid;
	theta_step  = options.theta_step;
	non_uniform = options.non_uniform;
	DISPLAY_EACH_ITERATION   = options.DISPLAY_EACH_ITERATION;
	SAVE_EACH_ITERATION   = options.SAVE_EACH_ITERATION;
	NUM_THREADS = options.NUM_THREADS;

	% Grid spacing for camera orientations
	theta_stepy = theta_step(1);
	theta_stepx = theta_step(2);
	theta_stepz = theta_step(3);

	% Camera orientations to which kernel elements correspond
	tty = theta_grid{1};
	ttx = theta_grid{2};
	ttz = theta_grid{3};

	Ksharp = Kblurry;

	imver = ver('images');
	if str2num(imver.Version) >= 5.4
		imresizefn = @imresize_old;
	else
		imresizefn = @imresize;
	end

	if non_uniform
		true_blur = reshape(true_blur,size(tty));
	end

	%%%%%%%%%%%%%
	%%%% Blur
	if strcmp(mode_blur,'direct')
		me2 = blur;
	elseif strcmp(mode_blur,'true')
		me2 = true_blur;
	elseif strcmp(mode_blur,'delta')
		[K,L] = size(true_blur);
		me2 = delta_kernel_uni(K);
	elseif strcmp(mode_blur,'hbar')
		[K,L] = size(true_blur);
		hK = floor(K/2)+1;
		hL = floor(L/2)+1;
		me2 = zeros(K,L);
		me2(hK,hL) = 1;
		me2(hK,hL-1) = 1;
		me2(hK,hL+1) = 1;
	elseif strcmp(mode_blur,'vbar')
		[K,L] = size(true_blur);
		hK = floor(K/2)+1;
		hL = floor(L/2)+1;
		me2 = zeros(K,L);
		me2(hK,hL) = 1;
		me2(hK-1,hL) = 1;
		me2(hK+1,hL) = 1;
	elseif strcmp(mode_blur,'star')
		[K,L] = size(true_blur);
		hK = floor(K/2)+1;
		hL = floor(L/2)+1;
		me2 = zeros(K,L);
		me2(hK-1,hL+1) = 1;
		me2(hK-1,hL-1) = 1;
		me2(hK+1,hL-1) = 1;
		me2(hK+1,hL+1) = 1;
	elseif strcmp(mode_blur,'random')
		[M,N] = size(true_blur);
		me2 = rand(M,N);
	elseif strcmp(mode_blur,'variational')
		me2 = blur;
		%%TODO
	else
		error foo
	end


	%%%%%%%%%%%%%%%
	%%%% Image

	spatial_mask = spatial_mask(:);

	if strcmp(mode_im,'direct')
		mx2 = im;
	elseif strcmp(mode_im,'true')
		mx2 = true_im;
	elseif strcmp(mode_im,'variational')
		MAX_ITERATIONS = 5000;
		[M,N] = size(im);
		[K,L] = size(blur);
		%norm_blur = blur / sum(blur(:));
		norm_blur = me2 / sum(me2(:));

		%%% n.b. for SIGGRAPH only 2 blur components were used.
		%%% set dimensions(2,6)=1 to get back to SSG & IMA runs

		dimensions = [1 1   1  0  0  1;
		              1 K*L 4  1  0  0;
		              1 M*N 4  0  1  1];



		I = M*2; J = N*2;  
		Dpf = zeros(I,J);
		if(~FFT_SHIFT)
			Dpf(K:M,L:N/2) = 1;
			Dpf(K:M,L+N/2:N) = 1; %% y-plane
		else
			hK = floor(K/2); hL = floor(L/2);
			Dpf(hK+1:M-hK,hL+1:N/2-hL) = 1;
			Dpf(hK+1:M-hK,N/2+hL+1:N-hL) = 1;
		end

		Df = padarray(obs,[M N],0,'post');

		pres_vector = ones(1,length(norm_blur(:))+length(im(:))+1) * pres; % precision
		q=find(spatial_mask(:)); % spatial_image_mask{s}
		pres_vector(q+1+length(norm_blur(:))) = spatial_mask(q)';

		dummy_blur_mask = zeros(dimensions(2,3),length(norm_blur(:)));

		%%% Make vectors
		xx1 = [0 norm_blur(:)' im(:)'] .* pres_vector;  
		xx2 = pres_vector;

		%     keyboard
		mainoptions = struct('converge_criteria',1e-4,...
							'plot_linesearch',0,...
							'DISPLAY_EACH_ITERATION',DISPLAY_EACH_ITERATION,...
							'SAVE_EACH_ITERATION',SAVE_EACH_ITERATION,...
							'init_noise_var',1,...
							'restart_priors',0,...
							'restart_switched_off',0,...
							'num_iters',MAX_ITERATIONS,...
							'text','',...
							'non_uniform',non_uniform,...
							'I',I,...
							'J',J,...
							'K',K,...
							'L',L,...
							'M',M,...
							'N',N,...
							'priors',priors,...
							'FFT_SHIFT',FFT_SHIFT,...
							'Kblurry',Kblurry,...
							'Ksharp',Ksharp,...
							'theta_grid',{theta_grid}, ...
							'NUM_THREADS',NUM_THREADS);
		[ensemble,D_log,gamma_log] = train_ensemble_main6(dimensions,xx1,xx2,mainoptions,Df,Dpf,dummy_blur_mask,[1-(spatial_mask>0)]);

		mx2 = reshape(train_ensemble_get(3,dimensions,ensemble.mx),M,N);

	end



	%%% ensure blur is normalized to 1
	summe2 = sum(me2(:));
	me2 = me2 / summe2;
	%  mx2  = mx2 * std(obs(:))/std(mx2(:));


	%%% do spatial masking
	pres_vector = ones(1,length(me2(:))+length(mx2(:))+nLayers) * pres;
	q=find(spatial_mask(:));
	pres_vector(q+1+length(me2(:))) = spatial_mask(q)';

	%%% Make vectors
	x1 = [zeros(1,nLayers) me2(:)' mx2(:)'] .* pres_vector;  
	x2 = pres_vector;
