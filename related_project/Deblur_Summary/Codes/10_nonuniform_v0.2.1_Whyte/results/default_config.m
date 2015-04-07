%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General configuration file, sets some defaults, parses some options
%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUTPUT_DIRECTORY = pwd;
imver = ver('images');
if str2num(imver.Version) >= 5.4
	imresizefn = @imresize_old;
else
	imresizefn = @imresize;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User specified parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin deblurring at a higher scale
if ~exist('START_SCALE','var')
	START_SCALE = 1;
end
% Default amount to downsample input image
if ~exist('PRESCALE','var')
	PRESCALE = 1;
end
% Actually run the kernel estimation, or load a previously saved kernel?
if ~exist('DO_ESTIMATE_KERNEL','var')
	DO_ESTIMATE_KERNEL = 1;
end
% Display a plot of the cost function along the search direction at each step
if ~exist('PLOT_LINESEARCH','var')
	PLOT_LINESEARCH = 0;
end
% Show the image & kernel estimated at each iterations
if ~exist('DISPLAY_EACH_ITERATION','var')
	DISPLAY_EACH_ITERATION = 0;
end
% Save the image & kernel estimated at each iterations
% WARNING: This will save thousands of images to disk
if ~exist('SAVE_EACH_ITERATION','var')
	SAVE_EACH_ITERATION = 0;
end
% Do deconvolution once kernel has been estimated?
if ~exist('DO_DECONVOLUTION','var')
	DO_DECONVOLUTION = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intrinsic calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('focal_length_35mm','var')
	focal_length_35mm = 35;
end

% downsample image before we start to keep image size managable
if (PRESCALE ~= 1)
    obs_im_scaled = imresizefn(obs_im,PRESCALE,'bilinear');
else
	obs_im_scaled = obs_im;
end

% Get size
[obs_imy,obs_imx,obs_imz] = size(obs_im_scaled);

% Calibration and range of angles
fl = max(size(obs_im_scaled))*focal_length_35mm/36;
x0 = obs_imx/2;
y0 = obs_imy/2;
Kblurry_all = [fl,0,x0;0,fl,y0;0,0,1];
max_radius = max(sqrt(([1, obs_imx, 1, obs_imx]-x0).^2 + ([1, 1, obs_imy, obs_imy]-y0).^2));

% Calibration and range of angles
flfr = max(size(obs_im))*focal_length_35mm/36;
x0fr = size(obs_im,2)/2;
y0fr = size(obs_im,1)/2;
Kblurry_fullres = [flfr,0,x0fr;0,flfr,y0fr;0,0,1];

clear obs_im_scaled

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Richardson - Lucy related parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of lucy richardson iterations to use
% default = 10 but more for long thin blurs...
if ~exist('LUCY_ITS','var')
	if NON_UNIFORM
		LUCY_ITS = 50;
	else
		LUCY_ITS = 50;
	end
end

% 0=use full res kernel, 1=use 1 scale coarser etc.
if ~exist('SCALE_OFFSET','var')
	SCALE_OFFSET = 0;
end

% threshold on kernel intensities (helps to remove noise from kernel)
% is % of max value in kernel
if ~exist('KERNEL_THRESHOLD','var')
	KERNEL_THRESHOLD = 7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other semi-important parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 = use rgb2gray; 1 = just use red channel, 2 = green; 3 = blue
if ~exist('SPECIFIC_COLOR_CHANNEL','var')
	SPECIFIC_COLOR_CHANNEL = 0; 
end

% type of prior to use
if ~exist('PRIOR_TYPE','var')
	% PRIOR_TYPE = 'whiteboard';
	PRIOR_TYPE = 'street';
end

% intensity value above which counts as saturation
if ~exist('SATURATION_THRESHOLD','var')
	SATURATION_THRESHOLD = 250;
end
% camera type 
if ~exist('CAMERA_TYPE','var')
	CAMERA_TYPE = 'unknown';
end

% gamma correction (set to 1 for RAW images)
if ~exist('GAMMA_CORRECTION','var')
	GAMMA_CORRECTION = 2.2;
end

% use automatic patch selector
if ~exist('AUTOMATIC_PATCH','var')
	AUTOMATIC_PATCH = 0;
end
% run synthetic experiment or not
if ~exist('SYNTHETIC','var')
	SYNTHETIC = 0;
end
% apply blur kernel priors or not
if ~exist('BLUR_LOCK','var')
	BLUR_LOCK = 1;
end
% re-center kernel after every scale
if ~exist('CENTER_BLUR','var')
	CENTER_BLUR = 1;
end

% Orientation quantization step
if ~exist('pixels_per_theta_step','var')
	pixels_per_theta_step = 1;
end

% for synthetic blur, synthesise it or load it from a file?
if SYNTHETIC && ~exist('EXTERNAL_SYNTHETIC','var')
	if NON_UNIFORM
		EXTERNAL_SYNTHETIC = 0;
	else
		EXTERNAL_SYNTHETIC = 1;
	end
end

if ~exist('KERNEL_DILATE_RADIUS','var')
	KERNEL_DILATE_RADIUS = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters - dont alter unless you really know what you are doing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RESIZE_STEP = sqrt(2); %% if ~=2 must use matlab..
NUM_SCALES = ceil(-log(3/BLUR_KERNEL_SIZE)/log(sqrt(2)))+1;

if ~AUTOMATIC_PATCH
  PATCH_SIZE = [AXIS(2)-AXIS(1) AXIS(4)-AXIS(3)];
  PATCH_LOCATION = [AXIS(1) AXIS(3)];
else
  PATCH_SIZE = [256 256];
end

BLUR_MASK = 0;
BLUR_MASK_VARIANCES = [50 0.02 5 0.2];

GRADIENT_MODE = 'haar';
%GRADIENT_MODE = 'steer';

%RESIZE_MODE = 'matlab_nearest';
RESIZE_MODE = 'matlab_bilinear';
%RESIZE_MODE = 'matlab_bicubic';
%RESIZE_MODE = 'binom5';

BLUR_MASK = 0;
BLUR_MASK_VARIANCES = [100 0.01];

AUTOMATIC_PATCH_CENTER_WEIGHT = 5; 
SATURATION_MASK = 1;

if (GAMMA_CORRECTION==1)
	warning('INTENSITY_SCALING = 1, change to 1/256 for 16-bit images?')
  INTENSITY_SCALING = 1;
else
  INTENSITY_SCALING = 1; 
end

if ~exist('RESCALE_THEN_GRAD','var')
RESCALE_THEN_GRAD = 0; %% order to do rescaling/gradient operations (only
                       %for non-eero functions
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUG = 0;
% SYNTHETIC = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization parameters to get x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRST_INIT_MODE_BLUR = 'lucy';
%FIRST_INIT_MODE_BLUR = 'reg';
%FIRST_INIT_MODE_BLUR = 'variational';
%FIRST_INIT_MODE_BLUR = 'direct';
%FIRST_INIT_MODE_BLUR = 'true';
%FIRST_INIT_MODE_BLUR = 'random';
%FIRST_INIT_MODE_BLUR = 'updown';
%FIRST_INIT_MODE_BLUR = 'vbar';
%FIRST_INIT_MODE_BLUR = 'delta';

%FIRST_INIT_MODE_IMAGE = 'lucy';
%FIRST_INIT_MODE_IMAGE = 'reg';
FIRST_INIT_MODE_IMAGE = 'variational';
%FIRST_INIT_MODE_IMAGE = 'direct';
%FIRST_INIT_MODE_IMAGE = 'true';
%FIRST_INIT_MODE_IMAGE = 'random';
%FIRST_INIT_MODE_IMAGE = 'slight_blur_obs';
%FIRST_INIT_MODE_IMAGE = 'updown';
%FIRST_INIT_MODE_IMAGE = 'nearest';
%FIRST_INIT_MODE_IMAGE = 'greenspan';

%INIT_MODE_BLUR = 'lucy';
%INIT_MODE_BLUR = 'reg';
%INIT_MODE_BLUR = 'variational';
INIT_MODE_BLUR = 'direct';
%INIT_MODE_BLUR = 'true';
%INIT_MODE_BLUR = 'updown';

%INIT_MODE_IMAGE = 'lucy';
%INIT_MODE_IMAGE = 'reg';
%INIT_MODE_IMAGE = 'variational';
INIT_MODE_IMAGE = 'direct';
%INIT_MODE_IMAGE = 'true';
%INIT_MODE_IMAGE = 'updown';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLUR_PRIOR = 1; % Exponential
IMAGE_PRIOR = 0; % Gaussian

if BLUR_LOCK
  BLUR_COMPONENTS = 4;
else
  BLUR_COMPONENTS = 3;
end

IMAGE_COMPONENTS = 4;
BLUR_UPDATE_FREQ = 1;
INIT_PRESCISION = 1e4;
if ~exist('NOISE_INIT')
	NOISE_INIT = 1;
end
if ~exist('CONVERGENCE')
	CONVERGENCE = 5e-4;
end
MAX_ITERATIONS = 50000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post iteration operations on blur kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%POST_PROCESS_BLUR = 'noise_removal';
POST_PROCESS_BLUR = 'none';
POST_PROCESS_THRESHOLD = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstsruct image mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('IMAGE_RECONSTRUCTION')
	IMAGE_RECONSTRUCTION = 'lucy';
end
%IMAGE_RECONSTRUCTION = 'reg';
%IMAGE_RECONSTRUCTION = 'variational';


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior file name root
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CAMERA_TYPE,'unknown')
  if strcmp(PRIOR_TYPE,'street')
    prior_name = ['../priors/linear_street_',num2str(IMAGE_COMPONENTS)];
  elseif strcmp(PRIOR_TYPE,'whiteboard')
    prior_name = ['../priors/linear_whiteboard_',num2str(IMAGE_COMPONENTS)];
  else
    error('Unknown prior type');
  end
else
  prior_name = ['../priors/',CAMERA_TYPE,'_street_',num2str(IMAGE_COMPONENTS)];
end
load(prior_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel for synthetic mode or blur guess for real mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SYNTHETIC && ~EXTERNAL_SYNTHETIC
	if NON_UNIFORM
		if ischar(SYNTH_ANGLE), SYNTH_ANGLE = eval(SYNTH_ANGLE); end
		theta_trajectory = shake_camera(zeros(obs_imy,obs_imx),Kblurry_all,[obs_imy,obs_imx],Kblurry_all,SYNTH_TYPE,SYNTH_ANGLE);
		blur_kernel = ones(size(theta_trajectory,2),1);
		blur_kernel = blur_kernel/sum(blur_kernel(:));
	else
		error('synthetic blur for uniform model must be loaded externally')
	end	
else
	if NON_UNIFORM
		blur_kernel = 'delta';
	else
		blur_kernel = delta_kernel_uni(BLUR_KERNEL_SIZE);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULT_GAMMA = 2.2;

if ~exist('MANUAL_KERNEL','var'),                         MANUAL_KERNEL = 0; end
if ~exist('FFT_SHIFT','var'),                                 FFT_SHIFT = 0; end
if ~exist('SATURATION_PRESCISION','var'),         SATURATION_PRESCISION = 1e1; end
if ~exist('CONSISTENCY_GRAD_ITS','var'),           CONSISTENCY_GRAD_ITS = 0; end
if ~exist('theta_trajectory','var'),                   theta_trajectory = []; end
if ~exist('THRESHOLD_SCALE_KERNEL','var'),       THRESHOLD_SCALE_KERNEL = 0; end % threshold kernel at end of each scale?
if ~exist('MASK_KERNEL','var'),                             MASK_KERNEL = 0; end % after thresholding, fix part of kernel to zero?
if ~exist('RESUME_PARTIAL','var'),                       RESUME_PARTIAL = 0; end % resume a partially completed result?
if ~exist('RESUME_PARTIAL_DEBLUR','var'),         RESUME_PARTIAL_DEBLUR = 0; end % resume a partially completed deblurring?
if ~exist('PRESCALE_BEFORE_GRAYSCALE','var'), PRESCALE_BEFORE_GRAYSCALE = 0; end % Fergus = 0
if ~exist('FERGUS_MASK_TRANSLATION','var'),		FERGUS_MASK_TRANSLATION = 1; end % Alternative mask translation to Fergus
% Number of threads for non-uniform model. Requires libpthread for NUM_THREADS > 1
if ~exist('NUM_THREADS','var'), 							NUM_THREADS = 1; end

% end of options section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save all to current filename
if SYNTHETIC
	OUTPUT_FILENAME = [CONFIG_FNAME,'_size',num2str(BLUR_KERNEL_SIZE),'_noise',num2str(SYNTHETIC_NOISE_STD)];
else
	OUTPUT_FILENAME = CONFIG_FNAME;
end
if ~exist(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME),'dir')
	mkdir(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME));
end

if DO_ESTIMATE_KERNEL
	save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'),'-v7.3');
else
	% Append, so we don't overwrite the file, it's got the kernel in it!
	save(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,'results.mat'),'-v7.3','-append');
end
