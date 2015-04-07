function out = delta_kernel_uni(siz)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% Modified by Oliver Whyte for CVPR 2010 paper:
% "Non-uniform Deblurring for Shaken Images"
% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

	%%% if even size, make odd
	if (mod(siz,2)==0)
		siz = siz + 1;
	end

	out = zeros(siz);

	c = floor(siz/2)+1;

	out(c,c)=1;


