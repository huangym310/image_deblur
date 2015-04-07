function [imout,pad] = edgetaper_uni(imblur,mask,psf)
% 	edgetaper_uni   blur edges of blurry image into zero with a non-uniform kernel
% 		[imout,pad] = edgetaper_uni(imblur,mask,psf)
% 	
% 	Author:		Oliver Whyte <oliver.whyte@ens.fr>
% 	Date:			August 2010
%	Copyright:		2010, Oliver Whyte
% 	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
% 	URL:			http://www.di.ens.fr/~whyte/deblurring/

% get imblur's dimensions
[h,w,c] = size(imblur);

% calculate maximum displacement
padx = size(psf,2);
pady = size(psf,1);

% make weights map
weights = double(mask);
% constraints = mask; constraints(:,[1,end]) = 1; constraints([1,end],:) = 1;
% weights = poisson_blend_backslash(zeros(size(mask)),zeros(size(mask)),mask,constraints);
weights = padarray(weights,[pady,padx],'both');
weights = repmat(weights,[1 1 size(imblur,3)]);

% pad blurred image
imblur = padarray(imblur,[pady,padx,0],'both');

fftpsf = fft2(ifftshiftpad(psf,[h,w]+2*[pady,padx]));

% blur imblur again
if c==3
	imtaper = real(cat(3,ifft2(fft2(imblur(:,:,1)).*fftpsf), ...
						 ifft2(fft2(imblur(:,:,2)).*fftpsf), ...
						 ifft2(fft2(imblur(:,:,3)).*fftpsf)));
else
	imtaper = real(ifft2(fft2(imblur).*fftpsf));
end
imtaper = max(imtaper,0);
% imtaper = apply_blur_kernel_mex(imblur,size(imblur(:,:,1))+2*[pady,padx],Kblurry,Kout,theta_list,blur_kernel/sum(blur_kernel(:)),0,1);

% crop to remove padding
pad = [pady,padx];

% combine images
imout = imtaper.*(1-weights) + imblur.*weights;

end %  function
