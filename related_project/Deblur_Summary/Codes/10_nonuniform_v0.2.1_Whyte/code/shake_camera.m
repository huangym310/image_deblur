function [theta_ss,im_blur] = shake_camera(im_sharp,Korig,dims_blur,Kblurry,blur_type,synth_size)
% 	shake_camera   Synthesize a camera-shake type blur
% 		[theta_ss,im_blur] = shake_camera(im_sharp,Korig,dims_blur,Kblurry,blur_type,synth_size)
% 
% 		Inputs:
% 				im_sharp   		sharp image
% 				Korig     		internal calibration matrix of sharp image
% 				dims_blur 		size of blurry image
% 				Kblurry     		internal calibration matrix of blurry image
% 				blur_type 		'random', 'x', 'y', or 'z'
% 				synth_size		size of blur to synthesize
% 
% 		Outputs:
% 				theta_ss		3 x 250 array containing trajectory of camera through orientations
% 									the 3 rows give coordinates for theta_X, theta_Y, theta_Z
% 				im_blur			blurry image
% 	
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

non_uniform = 1;

t = 0:0.01:0.1;
t_ss = linspace(min(t),max(t),250);
if strcmp(blur_type,'random')
    randn('seed',1);
    d2thetaxdt2 = randn(size(t)); d2thetaxdt2(1) = 0;
    d2thetaydt2 = randn(size(t)); d2thetaydt2(1) = 0;
    d2thetazdt2 = randn(size(t)); d2thetazdt2(1) = 0;
    theta0 = 1/100*[(cumsum(d2thetaxdt2)); ...
        			(cumsum(d2thetaydt2)); ...
        			(cumsum(d2thetazdt2))];
    theta_ss = interp1(t,theta0',t_ss,'spline')';
	theta_ss = theta_ss - repmat(mean(theta_ss,2),[1 size(theta_ss,2)]); % center kernel
	theta_ss = theta_ss/max(sqrt(sum(theta_ss.^2,1)))*synth_size*pi/180/2;
	theta_ss(3,:) = 2*theta_ss(3,:);
elseif strcmp(blur_type,'x')
    theta_ss = [linspace(-synth_size*pi/180/2,synth_size*pi/180/2,length(t_ss)); ...
		        zeros(size(t_ss)); ...
		        zeros(size(t_ss))];
elseif strcmp(blur_type,'y')
    theta_ss = [zeros(size(t_ss)); ...
		        linspace(-synth_size*pi/180/2,synth_size*pi/180/2,length(t_ss)); ...
		        zeros(size(t_ss))];
elseif strcmp(blur_type,'z')
    theta_ss = [zeros(size(t_ss)); ...
        		zeros(size(t_ss)); ...
        		linspace(-synth_size*pi/180/2,synth_size*pi/180/2,length(t_ss))];
end

if nargout > 1
	im_blur = zeros([dims_blur, size(im_sharp,3)]);
	for c=1:size(im_sharp,3)
		im_blur(:,:,c) = apply_blur_kernel_mex(im_sharp(:,:,c),dims_blur,Korig,Kblurry,-theta_ss,ones(size(theta_ss,2))/size(theta_ss,2),0,non_uniform,[0;0;0]);
	end
	im_blur = max(min(im_blur,1),0);
end


end %  function