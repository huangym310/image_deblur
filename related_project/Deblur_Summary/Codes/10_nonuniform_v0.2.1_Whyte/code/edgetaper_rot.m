function [imout,Kout,pad] = edgetaper_rot(imblur,mask,blur_kernel,theta_list,Kblurry)
% 	edgetaper_rot   Blur edges of blurry image into zero with a non-uniform kernel
% 		[imout,Kout,weights] = edgetaper_rot(imblur,mask,blur_kernel,theta_list,Kblurry)
% 
% 	Author:		Oliver Whyte <oliver.whyte@ens.fr>
% 	Date:			August 2010
%	Copyright:		2010, Oliver Whyte
% 	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
% 	URL:			http://www.di.ens.fr/~whyte/deblurring/

% get imblur's dimensions
[h,w,c] = size(imblur);

% calculate maximum displacement
corners = (theta_list(1,:)==max(theta_list(1,:)) | theta_list(1,:)==min(theta_list(1,:))) & ...
		  (theta_list(2,:)==max(theta_list(2,:)) | theta_list(2,:)==min(theta_list(2,:))) & ...
		  (theta_list(3,:)==max(theta_list(3,:)) | theta_list(3,:)==min(theta_list(3,:)));
corners = theta_list(:,corners);
imcorners = [1, 1, w, w; 1, h, 1, h; 1, 1, 1, 1];
disps = zeros(3,4*size(corners,2));
for cc=1:size(corners,2)
	disps(:,cc*4-3:cc*4) = hnormalise(Kblurry*expm(crossmatrix(corners(:,cc)))*inv(Kblurry)*imcorners) - imcorners;
end
padx = ceil(max(abs(disps(1,:))));
pady = ceil(max(abs(disps(2,:))));

% make output intrinsic calibration for padded image
Kout = htranslate([padx,pady])*Kblurry;

% make weights map
weights = double(mask);
% constraints = mask; constraints(:,[1,end]) = 1; constraints([1,end],:) = 1;
% weights = poisson_blend_backslash(zeros(size(mask)),zeros(size(mask)),mask,constraints);
weights = padarray(weights,[pady,padx],'both');
weights = repmat(weights,[1 1 size(imblur,3)]);

% blur imblur again
imtaper = apply_blur_kernel_mex(imblur,size(imblur(:,:,1))+2*[pady,padx],Kblurry,Kout,theta_list,blur_kernel/sum(blur_kernel(:)),0,1);

% pad blurred image
imblur = padarray(imblur,[pady,padx,0],'both');

% crop to remove padding
pad = [pady,padx];

% combine images
imout = imtaper.*(1-weights) + imblur.*weights;

end %  function