function out = delta_kernel_rot(siz)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% Modified by Oliver Whyte for CVPR 2010 paper:
% "Non-uniform Deblurring for Shaken Images"
% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

	out = zeros(siz);

	c = floor(prod(siz)/2)+1;

	out(c)=1;


