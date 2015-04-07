function L = non_blind_decov(B,K,lambda,iters,L)
% get image dimension
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
C = size(B,3);
eta = ones(C,1);
for c = 1:C
    eta_L = 1/2;
    for  it = 1:iters    
        [L(:,:,c),eta_L] = update_image_armijo(B(:,:,c),K,L(:,:,c),eta(c),lambda,eta_L);    
    end%it
end

end
