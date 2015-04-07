function plot_uni_kernel_write(K,file_name)
% 	plot_uni_kernel_write		Plot a non-uniform kernel and save it to disk
% 		[] = plot_uni_kernel_write(K,file_name)
% 
% 		Inputs:
% 				K				blur kernel
% 				file_name		file to write to
% 
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/


% Plot kernel
	h = figure('Visible','off'); 
	imagesc(K);
	colormap gray; axis off; axis image;
% Print to file
	print(h,'-dpng','-r100',file_name);
% Close figures
	close(h); 

