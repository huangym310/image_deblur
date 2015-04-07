function x = fftshiftpad(y,xdims)
% 	fftshiftpad   Zero-pad, then shift a filter y such that it has the same dimensions as the image to be convolved with it, and the central element is at (1,1)
% 		x = fftshiftpad(y,xdims)
% 
% 		Inputs:
% 				y		filter
% 				xdims	dimensions of image to be convolved with the filter
% 
% 		Outputs:
% 				x		shifted and zero-padded filter
% 	
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte

if nargin < 2, y = ifftshift(x); return; end

K = xdims(1);
L = xdims(2);

I = size(y,1);
J = size(y,2);

x = circshift(y,[ceil((K-1)/2), ceil((L-1)/2)]);
x = x(1:K,1:L);
