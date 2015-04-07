function T = htranslate(tt)
% 	htranslate   Translation matrix for homogeneous coordinates
% 		T = htranslate(tt)
% 
% 	For an input vector tt with n elements, the ouput T is the n+1 x n+1 matrix which applies that translation in homogeneous coordinates.
% 	
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte

T = eye(length(tt)+1);
T(1:end-1,end) = tt(:);

end %  function