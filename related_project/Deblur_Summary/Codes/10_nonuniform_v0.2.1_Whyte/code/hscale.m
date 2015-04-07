function S = hscale(ss)
% 	hscale   Scaling matrix for homogeneous coordinates
% 		S = hscale(ss)
% 
% 	For an input vector ss with n elements, S is a n+1 x n+1 matrix which applies those scaling factors in homogeneous coordinates.
% 	
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte


S = eye(length(ss)+1);
S(1:end-1,1:end-1) = diag(ss);

end %  function