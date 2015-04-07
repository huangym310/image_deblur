function K = kernel_normalize(K,threshold)
if ~exist('threshold','var')
    threshold = 0;
end
tau = max(K(:))*threshold;
K(K<tau) = 0;
K = K/sum(K(:));

end