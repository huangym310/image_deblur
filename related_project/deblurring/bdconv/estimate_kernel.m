function [kernel] = estimate_kernel(Ib,kernel_size,sample_rate)
%% input: blurred image Ib
%       parameters: 
%                  kernel_size -- 1-by-2 vector
%                  sample_rate -- sample_size = 2*sample_rate*kernel_size
%% output: blur kernel
%%
if nargin < 3 || isempty(sample_rate)
    sample_rate = 1.5;
end
[rsize_I,csize_I] = size(Ib);
if nargin<2 || isempty(kernel_size)
    rsize_ker = 2*round(rsize_I/50) + 1; 
    csize_ker = 2*round(csize_I/50) + 1;
else
    rsize_ker = kernel_size(1);
    csize_ker = kernel_size(2);
end
len_ker = rsize_ker*csize_ker;
disp('estimating the blur kernel ...');
%% 
%extract features 
rsize_sample = 2*round(sample_rate*rsize_ker) + 1;
csize_sample = 2*round(sample_rate*csize_ker) + 1;
g = fspecial('unsharp');
B = conv2(Ib,g,'same');
g = fspecial('log');
B = conv2(B,g,'same');
%% 
% compute the Hessian matrix
A = conv2mtx(B,rsize_sample,csize_sample,'same');
A = A'*A;
[U,S,~] = svd(A,'econ');
S = diag(S);
H = zeros(len_ker);
for i=1:size(U,2)
    Ai = U(:,i);
    Ai = reshape(Ai,rsize_sample,csize_sample);
    Ai = conv2mtx(Ai,rsize_ker,csize_ker,'same');
    H = H + (Ai'*Ai)/S(i);
end
%% 
% compute the kernel

f = zeros(len_ker,1);
options = optimset('quadprog');
options.Diagnostics = 'off';
options.Display = 'off';
options.Algorithm = 'active-set';
Aeq = ones(1, len_ker);
beq = 1;
lb = zeros(len_ker,1);
ub = ones(len_ker,1);
z = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
kernel = reshape(z,rsize_ker,csize_ker);
end






