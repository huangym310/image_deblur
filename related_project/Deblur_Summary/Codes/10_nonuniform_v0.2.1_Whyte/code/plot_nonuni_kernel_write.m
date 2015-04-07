function plot_nonuni_kernel_write(K,theta_grid,file_name,varargin)
% 	plot_nonuni_kernel_write		Plot a non-uniform kernel and save it to disk
% 		[] = plot_nonuni_kernel_write(K,theta_grid,file_name,varargin)
% 
% 		Inputs:
% 				K				blur kernel
% 				theta_grid		3 x 1 cell array, containing output of meshgrid
% 				file_name		file to write to
% 				varargin		other options, in pairs ... , 'name', [value] ...:
% 					 name			 value			 description
% 					========		========		==============
% 					'do_crop'		none			if present, crop kernel down to smallest cuboid enclosing non-zeros
% 					'view'			[az,el]			2 x 1 vector to set azimuth & elevation of camera
%
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

% Defaults
params.do_crop = 0;
% Get options
if nargin > 3, params = parseArgs(varargin,params); end
% Temporary files
tmpfile = tempname;
kernelfile = [tmpfile '_kernel.png'];
boxfile = [tmpfile '_kernelbox.png'];
% Actual kernel
h = figure('Visible','off'); plot_nonuni_kernel(K,theta_grid,1,0,1,params.do_crop);
if isfield(params,'view'); view(params.view(1),params.view(2)); end
set(h ,'InvertHardcopy','off')
print(h ,'-dpng','-r100',kernelfile);
% White box
h2 = figure('Visible','off'); plot_nonuni_kernel(K,theta_grid,1,1,1,params.do_crop);
if isfield(params,'view'); view(params.view(1),params.view(2)); end
set(h2,'InvertHardcopy','off')
print(h2,'-dpng','-r100',boxfile);
% Load white box and overlay onto kernel plot
I3 = max(imread(kernelfile), imread(boxfile));
% Extract alpha matte: red pixels are transparent
kernel_alpha = uint8(255*~(I3(:,:,1) == 255 & I3(:,:,2) == 0 & I3(:,:,3) == 0));
% Save image
imwrite(I3,file_name,'png','Alpha',kernel_alpha);
% Delete temporary files
delete(kernelfile);
delete(boxfile);
% Close figures
close(h); 
close(h2);
	
function out = parseArgs(args,out)
if ~isempty(args)
	switch args{1}
		case 'do_crop'
			out.do_crop = 1;
			args = args(2:end);
		case 'view'
			out.view = args{2};
			args = args(3:end);
		otherwise
			error('Invalid argument');
	end
	out = parseArgs(args,out);
end
