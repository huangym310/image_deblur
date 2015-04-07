function K = kernel_initialize(ks,kernel_type)
%output:
% K is ks X ks kernel
if ~exist('kernel_type','var')
    kernel_type = 'ver';
end
K = zeros(ks);
hks = floor(ks/2);
switch kernel_type
    case 'ver' %vertical
        K(hks(1)+1:hks(1)+2,hks(2)+1) = 1; 
    case 'hor' %horizontal
        K(hks(1)+1,hks(2)+1:hks(2)+2) = 1;
    case 'uniform'
        K(:) = 1;
    case 'delta'
        K(hks(1)+1,hks(2)+1) = 1;             
end% switch
K = K/sum(K(:));% normalize

end

