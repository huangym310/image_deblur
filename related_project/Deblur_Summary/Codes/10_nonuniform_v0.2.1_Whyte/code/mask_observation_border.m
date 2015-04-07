function mask = mask_observation_border(size_blurred,Kblurry,ttx,tty,ttz,Ksharp,size_sharp)
% 	mask_observation_border		Mask the edge of a blurry image, affected by boundary effects
% 		mask = mask_observation_border(size_blurred,Kblurry,ttx,tty,ttz,Ksharp,size_sharp)
% 
% 		Inputs:
% 				size_blurred		size of blurry image
% 				Kblurry       		internal calibration of blurry image
% 				ttx, tty, ttz		angles over which kernel varies, as produced by meshgrid
% 				Ksharp      		internal calibration of sharp image
% 				size_sharp  		size of sharp image
% 
% 		Outputs:
% 				mask				mask of pixels in blurry image unaffected by boundary effects
% 									1 = unaffected, 0 = affected by boundary
% 
%	Author:		Oliver Whyte <oliver.whyte@ens.fr>
%	Date:		August 2010
%	Copyright:	2010, Oliver Whyte
%	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
%	URL:		http://www.di.ens.fr/~whyte/deblurring/

    hei = size_sharp(1);
    wid = size_sharp(2);
    N = hei*wid;
    invKblur = inv(Kblurry);
    invKsharp = inv(Ksharp);
    heib = size_blurred(1);
    widb = size_blurred(2);

	% CORNERS OF CSF
    csfcorners = [[ttx(1,1,1), ttx(1,1,end), ttx(1,end,1), ttx(1,end,end), ttx(end,1,1), ttx(end,1,end), ttx(end,end,1), ttx(end,end,end)];...
                  [tty(1,1,1), tty(1,1,end), tty(1,end,1), tty(1,end,end), tty(end,1,1), tty(end,1,end), tty(end,end,1), tty(end,end,end)];...
                  [ttx(1,1,1), ttz(1,1,end), ttz(1,end,1), ttz(1,end,end), ttz(end,1,1), ttz(end,1,end), ttz(end,end,1), ttz(end,end,end)]];

	% FACES OF CSF
    % csfcorners = [[squash(ttx(1,:,:))', squash(ttx(end,:,:))', squash(ttx(:,1,:))', squash(ttx(:,end,:))', squash(ttx(:,:,1))', squash(ttx(:,:,end))'];...
    %               [squash(tty(1,:,:))', squash(tty(end,:,:))', squash(tty(:,1,:))', squash(tty(:,end,:))', squash(tty(:,:,1))', squash(tty(:,:,end))'];...
    %               [squash(ttx(1,:,:))', squash(ttz(end,:,:))', squash(ttz(:,1,:))', squash(ttz(:,end,:))', squash(ttz(:,:,1))', squash(ttz(:,:,end))']];

% p_test = [ x_test ; y_test ] is to the right of the oriented 
% line p_1 ---> p_2 if 0 < det([ (p_test-p_1) , (p_2-p_1) ])
%                        = (xt-x1)*(y2-y1) - (yt-y1)*(x2-x1)
%                        = xt*(y2-y1) + yt*(x1-x2) + x2*y1 - x1*y2
%                        = [xt;yt;1]' * [y2-y2;x1-x2;x2*y1-x1*y2]
%                        = [xt;yt;1]' * cross([x2;y2;1],[x1;y1;1])
mask = false(size_blurred);
imcorners2 = [1         ,1            ,size_sharp(2),size_sharp(2);...
              1         ,size_sharp(1),size_sharp(1),1            ;...
              1         ,1            ,1            ,1            ];
imedges2c1 = zeros(3,4*size(csfcorners,2));
for i=1:size(csfcorners,2)
         H = Kblurry*expm(crossmatrix(csfcorners(:,i)))*invKsharp;
         imcorners2c1 = inv(H)*imcorners2;
         imedges2c1(:,4*i-3:4*i) = [cross(imcorners2c1(:,1),imcorners2c1(:,4)), ... = [ y4-y1 ; x1-x4 ; x4*y1-x1*y4 ]
                       			    cross(imcorners2c1(:,2),imcorners2c1(:,1)), ...
                       				cross(imcorners2c1(:,3),imcorners2c1(:,2)), ...
                       				cross(imcorners2c1(:,4),imcorners2c1(:,3))];
end
% imedges2c1 = zeros(3,4*numel(ttx));
% for i=1:numel(ttx)
%          H = Kblurry*expm(crossmatrix([ttx(i);tty(i);ttz(i)]))*invKsharp;
%          imcorners2c1 = inv(H)*imcorners2;
%          imedges2c1(:,4*i-3:4*i) = [cross(imcorners2c1(:,1),imcorners2c1(:,4)), ... = [ y4-y1 ; x1-x4 ; x4*y1-x1*y4 ]
%                        				  cross(imcorners2c1(:,2),imcorners2c1(:,1)), ...
%                        				  cross(imcorners2c1(:,3),imcorners2c1(:,2)), ...
%                        				  cross(imcorners2c1(:,4),imcorners2c1(:,3))];
% end
[xx,yy] = meshgrid(1:size_blurred(2),1:size_blurred(1));    zz = ones(size_blurred);
mask = reshape(sum(sign([xx(:),yy(:),zz(:)] * imedges2c1),2) == 4*size(csfcorners,2),size_blurred);
    
% % ---------------------------------------------------------------------------------------
% % Calculate minimum necessary size of PSF to capture effect of shake at every pixel
% % ---------------------------------------------------------------------------------------
%     imcorners = [1  ,wid,1    ,wid;...
%                  1  ,1    ,hei,hei;...
%                  1  ,1    ,1    ,1    ];
%              todo fix this
% %     csfcorners = [[min(ttx(:)), min(ttx(:)), min(ttx(:)), min(ttx(:)), max(ttx(:)), max(ttx(:)), max(ttx(:)), max(ttx(:))];...
% %                   [min(tty(:)), min(tty(:)), min(tty(:)), min(tty(:)), max(tty(:)), max(tty(:)), max(tty(:)), max(tty(:))];...
% %                   [min(ttz(:)), min(ttz(:)), min(ttz(:)), min(ttz(:)), max(ttz(:)), max(ttz(:)), max(ttz(:)), max(ttz(:))]];
%     dxmin = inf; dxmax = -inf;
%     dymin = inf; dymax = -inf;
%     for i=1:size(csfcorners,2)
%         R = expm(crossmatrix(csfcorners(:,i)));
%         pp = hnormalise(Kblurry*R*invKsharp*imcorners);
%         dxmin = floor(min(dxmin,min(pp(1,:) - imcorners(1,:))));
%         dxmax =  ceil(max(dxmax,max(pp(1,:) - imcorners(1,:))));
%         dymin = floor(min(dymin,min(pp(2,:) - imcorners(2,:))));
%         dymax =  ceil(max(dymax,max(pp(2,:) - imcorners(2,:))));
%     end
%     dxmin = min(dxmin,0);    dxmax = max(dxmax,0);
%     dymin = min(dymin,0);    dymax = max(dymax,0);
% 
%     mask = zeros(size_blurred);
%     mask(dymax+1:heib+dymin,dxmax+1:widb+dxmin) = 1;
%     
%     
%     