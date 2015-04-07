function d = vec(image);
%
% Vectorize matrix "image"
% 
% d = vec(image)
%

i = prod(size(image));
d = reshape(image,i,1);
