function i_rl = deconvlucy_uni(imblur,psf,sat_thresh,num_iters)
% 	deconvlucy_uni		Apply Richardson-Lucy deconvolution with a uniform blur model
% 		i_rl = deconvlucy_uni(imblur,psf,sat_thresh,num_iters)
% 
% 		Inputs:
% 				imblur			Blurry image
% 				psf				Blur kernel to deconvolve with
% 				sat_thresh		Threshold above which, a pixel will be considered saturated
% 				num_iters		Number of Richarson-Lucy iterations
% 
% 		Outputs:
% 				i_rl			Deblurred image
% 
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

if nargin < 4, num_iters = 10; end
heio = size(imblur,1);
wido = size(imblur,2);
[heis,wids,channels] = size(imblur);

myscale = 1 + 255*(max(imblur(:))>1.1);

psf = psf/sum(psf(:));

% Calculate padding
hfs_x1 = 0;
hfs_x2 = 0;
hfs_y1 = 0;
hfs_y2 = 0;

hfs_x1 = ceil((size(psf,2)-1)/2);
hfs_x2 = floor((size(psf,2)-1)/2);
hfs_y1 = ceil((size(psf,1)-1)/2);
hfs_y2 = floor((size(psf,1)-1)/2);

wids=wids+hfs_x1+hfs_x2;
heis=heis+hfs_y1+hfs_y2;
wido=wido+hfs_x1+hfs_x2;
heio=heio+hfs_y1+hfs_y2;

% Replicate edges outwards
imblur = imblur([ones(1,hfs_y1),1:end,end*ones(1,hfs_y2)],[ones(1,hfs_x1),1:end,end*ones(1,hfs_x2)],:);

psf = ifftshiftpad(psf,[heis,wids]);
fftpsf = fft2(psf);
cftpsf = conj(fftpsf);

if channels == 3
	blurfn = @(im) real(cat(3,ifft2(fft2(im(:,:,1)).*fftpsf),ifft2(fft2(im(:,:,2)).*fftpsf),ifft2(fft2(im(:,:,3)).*fftpsf)));
	blurconjfn = @(im) real(cat(3,ifft2(fft2(im(:,:,1)).*cftpsf),ifft2(fft2(im(:,:,2)).*cftpsf),ifft2(fft2(im(:,:,3)).*cftpsf)));
else
	blurfn = @(im) real(ifft2(fft2(im).*fftpsf));
	blurconjfn = @(im) real(ifft2(fft2(im).*cftpsf));
end

imblur = imblur/myscale;

% Initialise sharp image by projecting blurry image into sharp frame
i_rl = imblur;

for iter = 1:num_iters
    i_rl = blurconjfn(imblur./max(blurfn(i_rl),eps)) .* i_rl;
end

i_rl = i_rl(hfs_y1+1:end-hfs_y2,hfs_x1+1:end-hfs_x2,:)*myscale;


