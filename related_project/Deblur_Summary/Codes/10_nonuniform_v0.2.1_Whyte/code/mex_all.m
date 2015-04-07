% Mex all the functions in the deblurring code

mex apply_blur_kernel_mex.c
mex mask_observation_border_mex.c
mex train_blind_deconv_mex.c

% Requires libpthread, so far Unix only
if isunix
	mex train_blind_deconv_mex_threads.c -lpthread
end
