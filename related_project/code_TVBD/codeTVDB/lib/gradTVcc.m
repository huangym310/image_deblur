function divTV = gradTVcc(f,epsilon)
% Total Variation gradient
% vector-valued channel-by-channel
%
% Author: Paolo Favaro
% Copyright (C) 2013 Paolo Favaro. All rights reserved.

fxforw = f([2:end end],:,:)-f;
fyforw = f(:,[2:end end],:)-f;
fxback = f-f([1 1:end-1],:,:);
fyback = f-f(:,[1 1:end-1],:);
fxmixd = f([2:end end],[1 1:end-1],:)-f(:,[1 1:end-1],:);
fymixd = f([1 1:end-1],[2:end end],:)-f([1 1:end-1],:,:);
if nargin<2
    epsilon = 1e-3;
end
divTV = zeros(size(f));
for cc=1:size(f,3)
    divTV(:,:,cc)  = ...
        ( fxforw(:,:,cc) + fyforw(:,:,cc))./max(epsilon,sqrt(fxforw(:,:,cc).^2+fyforw(:,:,cc).^2))...
         -fxback(:,:,cc)                  ./max(epsilon,sqrt(fxback(:,:,cc).^2+fymixd(:,:,cc).^2))...
                         - fyback(:,:,cc) ./max(epsilon,sqrt(fxmixd(:,:,cc).^2+fyback(:,:,cc).^2));
end
