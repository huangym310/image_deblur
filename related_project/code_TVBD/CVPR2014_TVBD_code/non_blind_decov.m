function L = non_blind_decov(B,K,lambda,iters,L)
% get image dimension
[M,N,C] = size(B);
[MK,NK] = size(K);
B = im2double(B);
if ~exist('lambda','var')
    lambda = 3e-4;
end
if ~exist('L','var')
    L = padarray(B,[floor(MK/2),floor(NK/2)],'replicate');
end
if ~exist(' iters','var')
    iters = 1000;
end

% size of sharp image
ML = M + MK-1;
NL = N + NK-1;
gradudata = zeros(ML,NL,C);
for  it = 1:iters
    for c = 1:C        
        gradudata(:,:,c) = conv2(conv2(L(:,:,c),K,'valid')-B(:,:,c),...
                                 rot90(K,2), 'full');
    end   
    dL = gradudata+lambda*gradTVcc(L);
    eta_L = 5e-3*max(L(:))/max(1e-31,max(abs(dL(:))));
    L = L-eta_L*dL;    
end%it

end
