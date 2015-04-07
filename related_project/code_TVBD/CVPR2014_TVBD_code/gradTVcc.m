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
divTV  = -(fxforw+fyforw)./max(epsilon,sqrt(fxforw.^2+fyforw.^2))...
         +fxback./max(epsilon,sqrt(fxback.^2+fymixd.^2))...
         +fyback./max(epsilon,sqrt(fxmixd.^2+fyback.^2));

end%function
