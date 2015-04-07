function vx = crossmatrix(v)
% 	crossmatrix		form 3 x 3 matrix equivalent to vector cross product
%		vx = crossmatrix(v)

%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

vx = [    0, -v(3),  v(2);...
       v(3),     0, -v(1);...
      -v(2),  v(1),     0];