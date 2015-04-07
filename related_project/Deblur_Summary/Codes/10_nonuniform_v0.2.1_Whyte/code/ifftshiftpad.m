function y = ifftshiftpad(x,ydims)
% 	ifftshiftpad   Undo the shift and zero-padding of fftshiftpad
% 		y = ifftshiftpad(x,ydims)
% 
% 		Inputs:
% 				x		filter, zero-padded and shifted to have its central element at (1,1)
% 				ydims	dimensions of the original filter
% 
% 		Outputs:
% 				y		de-padded and unshifted filter
% 	
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte

if nargin < 2, y = ifftshift(x); return; end

K = size(x,1);
L = size(x,2);

I = ydims(1);
J = ydims(2);

y = padarray(x,[I-K,J-L],'post');
y = circshift(y,[-ceil((K-1)/2), -ceil((L-1)/2)]);
