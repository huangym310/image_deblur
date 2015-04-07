SYNTHETIC = 1;
SYNTH_TYPE = 'z';

if NON_UNIFORM
	blur_x_lims = floor((sqrt(BLUR_KERNEL_SIZE)-1)/2)*[-1 1];
	blur_y_lims = floor((sqrt(BLUR_KERNEL_SIZE)-1)/2)*[-1 1];
	blur_z_lims = floor((BLUR_KERNEL_SIZE-1)/2)*[-1 1];
	% inital value of kernel
	FIRST_INIT_MODE_BLUR = 'zbar';
	SYNTH_ANGLE = '(blur_z_lims(2)-blur_z_lims(1))/max_radius*0.75*180/pi';
else
	EXTERNAL_SYNTHETIC = 1;
	% inital value of kernel
	FIRST_INIT_MODE_BLUR = 'hbar';
end

GAMMA_CORRECTION = 1;
