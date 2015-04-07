function [var1 var2 var3] = homography(H,igridX,igridY,I,varargin)
%
% transform mesh igrid using homography H
%
% homography(H,igridX,igridY,I,method,extrapval)
%
% igrid ... input mesh, if empty default mesh is used
% H ... homography 3x3 matrix
% I ... image on which we apply H
% method ... interpolation method (see interp2)
% extrapval ... extrapolation value (see interp2)
%
% variable number of outputs
% ogrid ... output grid
% R ... transformed image
%

% if input grid is empty use the default one
if isempty(igridX)
   [igridX, igridY] = meshgrid((1:size(I,2))-ceil(size(I,2)/2),(1:size(I,1))-ceil(size(I,1)/2));
end
A = H*[igridX(:),igridY(:),ones(numel(igridX(:)),1)].';
ogridX = reshape(A(1,:)./A(3,:),size(igridX));
ogridY = reshape(A(2,:)./A(3,:),size(igridY));

if nargin >= 4
   R = zeros(size(ogridX,1),size(ogridX,2),size(I,3));
   for i = 1:size(I,3) 
      R(:,:,i) = interp2((1:size(I,2))-ceil(size(I,2)/2),(1:size(I,1))-ceil(size(I,1)/2),...
        I(:,:,i),ogridX,ogridY,varargin{:});
   end
   if nargout == 1
      var1 = R;
   else
      var1 = ogridX;
      var2 = ogridY;
      var3 = R;
   end
else
   var1 = ogridX;
   var2 = ogridY;
end
    