function plot_nonuni_kernel(K,theta_grid,currfig,show_box,plotpoints,do_crop,friendly)
% 	plot_nonuni_kernel	Function to plot camera-rotation kernels, as in cvpr10 paper.
% 		[] = plot_nonuni_kernel(K,theta_grid,currfig,show_box,plotpoints,do_crop,friendly)
% 	
% 		Inputs:
% 				K				blur kernel
% 				theta_grid		3 x 1 cell array of angles K covers, containing output of meshgrid
% 				currfig			if true, plot in current figure, else make new figure
% 				show_box		if true, draw white bounding box, else plot the kernel
% 				plotpoints		if true, plot non-zeros as yellow points in 3D
% 				do_crop			if true, crop kernel down to smallest cuboid enclosing the non-zeros
% 				friendly		if true, don't set any figure properties, e.g. font size / background colour
% 
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/


if nargin < 3 || isempty(currfig), currfig = 0; end
if nargin < 4 || isempty(show_box), show_box = 0; end
if ~currfig
	figure;
end
if nargin < 5 || isempty(plotpoints), plotpoints = 0; end
if nargin < 6 || isempty(do_crop), do_crop = 0; end
% 'friendly' if we just want to visualize the kernel during an experiment, not necessarily save it out to a file. So make it look nice and don't do the crazy background colour changes, or set the font to be giant.
if nargin < 7 || isempty(friendly), friendly = 0; end

force_3d = 0;

cla; hold off

tgs = median(diff(sort(unique(theta_grid{1}))));
if isnan(tgs)
	tgs = median(diff(sort(unique(theta_grid{2}))));
end
tgsz = median(diff(sort(unique(theta_grid{3}))));
if isnan(tgsz)
	tgsz = tgs/2;
end

Kpn = reshape(K,size(theta_grid{1}));

% Expand any singular dimensions
ksize = ones(1,3);
ksize(1:ndims(Kpn)) = size(Kpn);
sdims = ksize == 1;

expfully = 1; % expand to be a full cuboid? if not, just expand singular dimensions by 1 on each side
% maxexp = max(ksize);
maxexp = 3; % if expfully, maxep = minumum allowable dimension along each axis. below this, a dimension will be expanded to maxexp 

% which dimensions to expand?
if expfully
	expdims = find(ksize < maxexp);
else
	expdims = find(sdims);
end
% if we don't want to go to 3d when kernel is 1d or 2d, remove last entry of expdims
if any(sdims) && ~force_3d
	dontexpdim = expdims(end); % which dimension we don't expand
	expdims(end) = [];
end

tt0 = theta_grid;

for sdim = expdims
	% Extend theta arrays:
	repmatdims = [1, 1, 1];
	if expfully, repmatdims(sdim) = maxexp; else, repmatdims(sdim) = 3; end
	halfrepmatdims = repmatdims;    halfrepmatdims(sdim) = (repmatdims(sdim)-size(Kpn,sdim))/2;
	% Pad kernel along singular dimension with a layer of zeros either side
	zerodims = size(Kpn); zerodims(sdim) = halfrepmatdims(sdim);
	Kpn = cat(sdim,zeros(zerodims),Kpn,zeros(zerodims));
	% For singular dimension, extrapolate values
	tgsadd = reshape((1:halfrepmatdims(sdim)),halfrepmatdims);
	if sdim == 3
		tgsadd = tgsz*tgsadd;
		ttslice1 = theta_grid{sdim}(:,:,1);
		ttslice2 = theta_grid{sdim}(:,:,end);
		% theta_grid{sdim} = cat(sdim,theta_grid{sdim}-tgsz,theta_grid{sdim},theta_grid{sdim}+tgsz);
	elseif sdim == 2
		tgsadd = tgs*tgsadd;
		ttslice1 = theta_grid{sdim}(:,1,:);
		ttslice2 = theta_grid{sdim}(:,end,:);
		% theta_grid{sdim} = cat(sdim,theta_grid{sdim}-tgs, theta_grid{sdim},theta_grid{sdim}+tgs );
	else
		tgsadd = tgs*tgsadd;
		ttslice1 = theta_grid{sdim}(1,:,:);
		ttslice2 = theta_grid{sdim}(end,:,:);
	end
	% theta_grid{sdim} = cat(sdim, bsxfun(@minus, repmat(ttslice1, halfrepmatdims), tgsadd(end:-1:1)), theta_grid{sdim}, bsxfun(@plus, repmat(ttslice2, halfrepmatdims), tgsadd));
	ttslice1exp = repmat(ttslice1, halfrepmatdims);
	ttslice2exp = repmat(ttslice2, halfrepmatdims);
	tgssubexp = repmat(tgsadd(end:-1:1),[size(ttslice1exp,1),size(ttslice1exp,2),size(ttslice1exp,3)]./[size(tgsadd,1),size(tgsadd,2),size(tgsadd,3)]);
	tgsaddexp = repmat(tgsadd          ,[size(ttslice1exp,1),size(ttslice1exp,2),size(ttslice1exp,3)]./[size(tgsadd,1),size(tgsadd,2),size(tgsadd,3)]);
	% keyboard
	theta_grid{sdim} = cat(sdim, ttslice1exp-tgssubexp, theta_grid{sdim}, ttslice2exp+tgsaddexp);
	% For other dimensions, replicate values
	for nsdim = find([1, 2, 3] ~= sdim)
		if sdim == 1
			ttslice = theta_grid{nsdim}(1,:,:);
		elseif sdim == 2
			ttslice = theta_grid{nsdim}(:,1,:);
		else
			ttslice = theta_grid{nsdim}(:,:,1);
		end
		theta_grid{nsdim} = repmat(ttslice,repmatdims);
		% theta_grid{nsdim} = cat(sdim,theta_grid{nsdim},)
	end
end	
ttxn = theta_grid{2}*180/pi;
ttyn = theta_grid{1}*180/pi;
ttzn = theta_grid{3}*180/pi;

% get kernel dimensions now
ksize = ones(1,3);
ksize(1:ndims(Kpn)) = size(Kpn);
sdims = ksize == 1;
nkdims = nnz(~sdims);

Kpn = Kpn/max(Kpn(:));

kpy = sum(Kpn,1);   kpy = squeeze(kpy)/max(kpy(:));
kpx = sum(Kpn,2);   kpx = squeeze(kpx)/max(kpx(:));
kpz = sum(Kpn,3);   kpz = squeeze(kpz)/max(kpz(:));

% Crop kernel to center non-zeros?
if do_crop
	cropmargin = 5;
	[nzy,nzx,nzz] = ind2sub(ksize,find(Kpn~=0));
	centre_y = round(mean(nzy(:)));
	centre_x = round(mean(nzx(:)));
	centre_z = round(mean(nzz(:)));
	range_y = [max(1,min(nzy(:))-cropmargin), min(ksize(1),max(nzy(:))+cropmargin)];
	range_x = [max(1,min(nzx(:))-cropmargin), min(ksize(2),max(nzx(:))+cropmargin)];
	range_z = [max(1,min(nzz(:))-cropmargin), min(ksize(3),max(nzz(:))+cropmargin)];
	% Make range_x and range_y match as best as possible
	range_y = [max(1,min(range_y(1), centre_y-(centre_x-range_x(1)))), min(ksize(1),max(range_y(2), centre_y+(range_x(2)-centre_x)))];
	range_x = [max(1,min(range_x(1), centre_x-(centre_y-range_y(1)))), min(ksize(2),max(range_x(2), centre_x+(range_y(2)-centre_y)))];
	% Make limits into ranges
	range_y = range_y(1):range_y(2);
	range_x = range_x(1):range_x(2);
	range_z = range_z(1):range_z(2);
	% Do cropping
	Kpn = Kpn(range_y,range_x,range_z);
	ttxn = ttxn(range_y,range_x,range_z);
	ttyn = ttyn(range_y,range_x,range_z);
	ttzn = ttzn(range_y,range_x,range_z);
	kpy = kpy(range_x,range_z);
	kpx = kpx(range_y,range_z);
	kpz = kpz(range_y,range_x);
end

% color issues
m = 256;
cm = gray(m);
% cm = cm.^4;
% cm = cm(end:-1:1,:); % reverse colormap -> white to black increasing
cminx = min(kpx(:));   cmaxx = max(kpx(:));   idxx = min(m,round((m-1)*(kpx-cminx)/(cmaxx-cminx))+1);
cminy = min(kpy(:));   cmaxy = max(kpy(:));   idxy = min(m,round((m-1)*(kpy-cminy)/(cmaxy-cminy))+1);
cminz = min(kpz(:));   cmaxz = max(kpz(:));   idxz = min(m,round((m-1)*(kpz-cminz)/(cmaxz-cminz))+1);

hold on

if nkdims == 3
	hh1 = warp(squeeze(ttxn(1,:,:)),  squeeze(ttyn(1,:,:)),  squeeze(ttzn(1,:,:)),  reshape(cm(idxy,:),[size(kpy) 3]));
	hh2 = warp(squeeze(ttxn(:,end,:)),squeeze(ttyn(:,end,:)),squeeze(ttzn(:,end,:)),reshape(cm(idxx,:),[size(kpx) 3]));
	hh3 = warp(squeeze(ttxn(:,:,1)),  squeeze(ttyn(:,:,1)),  squeeze(ttzn(:,:,1)),  reshape(cm(idxz,:),[size(kpz) 3]));

	nz_ttx = ttxn(Kpn(:)>0); nz_tty = ttyn(Kpn(:)>0); nz_ttz = ttzn(Kpn(:)>0);

	CData = Kpn(Kpn(:)>0);
	cmax = max(CData(:));
	cmin = min(CData(:));
	if cmin==cmax, cmin=0; end
	idx = min(m,round((m-1)*(CData-cmin)/(cmax-cmin))+1);
	% scatter3(nz_ttx,nz_tty,nz_ttz,Kpn(Kpn(:)>0)*0+50,cm(idx,:),'filled');
	if plotpoints
		hh4 = scatter3(nz_ttx,nz_tty,nz_ttz,Kpn(Kpn(:)>0)*0+100,'y','filled');
		hh5 = scatter3(nz_ttx,nz_tty,nz_ttz,Kpn(Kpn(:)>0)*0+100,'k','LineWidth',1.5);
	end
	axis square;
	axis tight;
	% axis off;
	% currax=axis; axx=[-1 1]*max(abs(currax(1:2))); axy=[-1 1]*max(abs(currax(3:4))); axz=[-1 1]*max(abs(currax(5:6))); axis([axx axy axz])
	%currax=axis; axx=[-1 1]*max(abs(ttxn(:))); axy=[-1 1]*max(abs(ttyn(:))); axz=[-1 1]*max(abs(ttzn(:))); axis([axx axy axz])
	camproj('perspective'); % set(gca,'Color',[0 0 0]);
	if show_box
		% delete all plotted things
		v = axis;
		delete(hh1);
		delete(hh2);
		delete(hh3);
		set(gcf,'renderer','opengl'); % otherwise reverts to painters once the warps are deleted
		axis(v);
		% delete(hh4);
		% delete(hh5);
		% hh6 = warp(squeeze(ttxn(end,:,:)),  squeeze(ttyn(end,:,:)),  squeeze(ttzn(end,:,:)),  reshape(cm(idxy*0+1,:),[size(kpy) 3]));
		% hh7 = warp(squeeze(ttxn(:,1,:)),squeeze(ttyn(:,1,:)),squeeze(ttzn(:,1,:)),reshape(cm(idxx*0+1,:),[size(kpx) 3]));
		% hh8 = warp(squeeze(ttxn(:,:,end)),  squeeze(ttyn(:,:,end)),  squeeze(ttzn(:,:,end)),  reshape(cm(idxz*0+1,:),[size(kpz) 3]));
	end
	% aspect ratio of plot, need to swap x and y
	myaspect = size(Kpn);   myaspect([1,2]) = myaspect([2,1]);
	pbaspect(myaspect); % set aspect of plot-box so that voxels are cubes
else
	if sdims(1)
		imagesc(reshape(cm(idxy,:),[size(kpy) 3]));
	elseif sdims(2)
		imagesc(reshape(cm(idxx,:),[size(kpx) 3]));
	elseif sdims(3)
		% rotate so theta_y is horizontal
		imagesc(permute(reshape(cm(idxz,:),[size(kpz) 3]),[2,1,3]));
	end
	axis square;
	axis tight;
	camproj('orthographic'); % set(gca,'Color',[0 0 0]);
	view(0,90);
end
% show axis above projections
set(gca,'Layer','top','XTick',[],'YTick',[],'ZTick',[]);

if show_box
	set(gca,'Color','none','Box','on','LineWidth',3,'XColor','white','YColor','white','ZColor','white');
	set(gcf,'Color','black');
else
	if ~friendly
		set(gcf,'Color','red');
	end
end

if nkdims == 3
	% text axis labels
	xlabel('\theta_X','Color','black');
	ylabel('\theta_Y','Color','black');
	zlabel('\theta_Z','Color','black')
else
	if sdims(1)
		xlabel('\theta_X','Color','black');
		ylabel('\theta_Z','Color','black');
	elseif sdims(2)
		xlabel('\theta_Y','Color','black');
		ylabel('\theta_Z','Color','black');
	elseif sdims(3)
		xlabel('\theta_Y','Color','black');
		ylabel('\theta_X','Color','black');
	end
end
if ~friendly
	set_figure_font_size(gcf,40);
end

if expfully && nkdims == 3
	view(-40,40); % set viewpoint so lines don't coincide
end

drawnow
