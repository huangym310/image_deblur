% User selected region (xmin xmax ymin ymax)
AXIS = [1 259 1 194];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image to deblur & intrinsic calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename of image to deblur
obs_im = imread('../images/pantheon.jpg');

% Focal length in 35mm equivalent
focal_length_35mm = 35.1;

% downsample image before we start?
PRESCALE = 0.1;

% size of blur kernel in the image
BLUR_KERNEL_SIZE = 9;

% inital value of kernel
% FIRST_INIT_MODE_BLUR = 'hbar';
FIRST_INIT_MODE_BLUR = 'vbar';
%FIRST_INIT_MODE_BLUR = 'delta';

% parameters for dimensions of non-uniform kernel
pixels_per_theta_step = 1;
blur_x_lims = floor(((BLUR_KERNEL_SIZE)-1)/2)*[-1 1];
blur_y_lims = floor(((BLUR_KERNEL_SIZE)-1)/2)*[-1 1];
blur_z_lims = floor(((BLUR_KERNEL_SIZE)-1)/4)*[-1 1];


