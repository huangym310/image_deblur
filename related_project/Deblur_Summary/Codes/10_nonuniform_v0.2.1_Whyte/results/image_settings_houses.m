% User selected region (xmin xmax ymin ymax)
AXIS = [1 320 1 240];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image to deblur & intrinsic calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename of image to deblur
obs_im = imread('../images/houses.jpg');  

% downsample image before we start to keep image size managable
PRESCALE = 1;

% Focal length in 35mm equivalent
focal_length_35mm = 36;

SATURATION_THRESHOLD = 1000;

GAMMA_CORRECTION = 1;

DO_ESTIMATE_KERNEL = 1;
