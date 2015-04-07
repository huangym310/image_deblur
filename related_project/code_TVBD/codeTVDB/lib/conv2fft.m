function c = conv2fft(a,b,domain)
% Author: Paolo Favaro, paolo.favaro@iam.unibe.ch
% Copyright (C) 2013 Paolo Favaro.  All rights reserved.

[Nx1,Nx2] = size(a);
[NKx1,NKx2] = size(b);
ahat = fftshift(fftn(a,[Nx1+NKx1-1 Nx2+NKx2-1]));
bhat = fftshift(fftn(b,[Nx1+NKx1-1 Nx2+NKx2-1]));

chat = ahat.*bhat;
% full domain
c = real(ifftn(ifftshift(chat)));

if strcmp(domain,'same')
    c = c(floor(NKx1/2)+[1:Nx1],floor(NKx2/2)+[1:Nx2]);
elseif strcmp(domain,'valid')
    c = c(NKx1-1+[1:Nx1-NKx1+1],NKx2-1+[1:Nx2-NKx2+1]);
end
return