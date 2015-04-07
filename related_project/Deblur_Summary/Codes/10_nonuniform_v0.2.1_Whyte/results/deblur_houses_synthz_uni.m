if ~exist('BLUR_KERNEL_SIZE','var') || ~exist('SYNTHETIC_NOISE_STD','var')
	error('Must set BLUR_KERNEL_SIZE and SYNTHETIC_NOISE_STD before running this script')
end

CONFIG_FNAME = 'houses_synthz_uni';

NON_UNIFORM = 0;

image_settings_houses;

synthetic_settings_theta_z;

default_config;

deblur;