function [D,rankD] = sv_shrink(W,epsilon)
% single value shrinkage
% solve the following optimization problem:
% D = arg_X  min epsilon||X||_* + 1/2||X-W||_F^2

[U,S,V] = svd(W,'econ');%S is a square matrix, the size is min(m,n)            
Svec = diag(S);
S_epsilon = shrink(Svec,epsilon);%soft-thresholding
idx = (S_epsilon~=0);
D = U(:,idx)*diag(S_epsilon(idx))*V(:,idx)';
if nargout > 1
    rankD = sum(Svec>epsilon);%the rank of D
end

end%of function