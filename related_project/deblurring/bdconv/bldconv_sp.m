function [kernel,I0] = bldconv_sp(Ib,kernel_size,sample_rate,alpha,lambda)
%% input: blurred image Ib (we assume that the pixel values have been normlized to [0,1])
%       parameters: 
%                  kernel_size -- 1-by-2 vector
%                  sample_rate -- sample_size = 2*sample_rate*kernel_size
%                  alpha -- regularizer for kernel prior
%                  lambda -- regularizer for image prior
%% output:  blur kernel, deblured image
%%
[rsize_I,csize_I] = size(Ib);
rsize_ker = kernel_size(1);
csize_ker = kernel_size(2);
len_ker = rsize_ker*csize_ker;
rsize_sample = 2*round(sample_rate*rsize_ker) + 1;
csize_sample = 2*round(sample_rate*csize_ker) + 1;
rrad_ker = round((rsize_ker-1)/2);
crad_ker = round((csize_ker-1)/2);
saveres = false;
if saveres
    sdir = ['temp/' num2str(round(rand()*100000000))];
    mkdir(sdir);
end
%% 
disp('computing the Hessian matrix ...');
%extract features 
g = fspecial('unsharp');
B = conv2(Ib,g,'same');
g = fspecial('log');
B = conv2(B,g,'same');
%% 
A = conv2mtx(B,rsize_sample,csize_sample,'same');
% compute the Hessian matrix
A = A'*A;
if issparse(A)
    A = full(A);
end
[U,S,V] = svd(A,'econ');
S = diag(S);
H = zeros(len_ker);
for i=1:size(U,2)
    Ai = U(:,i);
    Ai = reshape(Ai,rsize_sample,csize_sample);
    Ai = conv2mtx(Ai,rsize_ker,csize_ker,'same');
    H = H + (Ai'*Ai)/S(i);
end
%
%% loop
%     
I0 = Ib;

rm1 = min(2*rsize_ker,round(rsize_I/5));
cm1 = min(2*csize_ker,round(csize_I/5));
rm = rm1+rrad_ker;
cm = cm1+crad_ker;
vB = Ib(rm + 1:rsize_I - rm, cm + 1:csize_I - cm);
vB = reshape(vB,size(vB,1)*size(vB,2),1);

options = optimset('quadprog');
options.Diagnostics = 'off';
options.Display = 'off';
options.LargeScale = 'off';
options.Algorithm = 'active-set';
Aeq = ones(1,len_ker);
beq = 1;

lb = zeros(len_ker,1);
ub = ones(len_ker,1);
x0 = ub./sum(ub);
disp('loop ...');
iter = 0;
lambda = 1/lambda;
for level=1:6
    if level == 1
        max_iter = 10;
        lambda1 = lambda/20;
    elseif level == 2
        max_iter = 20;
        lambda1 = lambda/15;
    elseif level ==3 
        max_iter = 35;
        lambda1 = lambda/10;
    elseif level == 4
        max_iter = 50;
        lambda1 = lambda/7;
    elseif level ==5 
        max_iter = 70;
        lambda1 = lambda/3;
    else
        max_iter = 200;
        lambda1 = lambda;
    end
    obj0 = inf;
    while iter < max_iter
        iter = iter + 1;

        I01 = I0(rm1+1:rsize_I-rm1, cm1+1:csize_I-cm1);
        A = conv2mtx(I01,rsize_ker,csize_ker,'valid');
        f =  - A'*vB;
        G = alpha*H + A'*A;

        [z,obj] = quadprog(G,f,[],[],Aeq,beq,lb,ub,x0,options);

        kernel = reshape(z,rsize_ker,csize_ker);
        
        I0 = fast_deconv(Ib, kernel, lambda1, 1, I0);
        
        obj = obj + 0.5*norm(vB).^2;

        disp([' iter ' num2str(iter) ', obj=' num2str(obj)]);
        if abs(obj - obj0) < 1e-6*obj 
            break;
        else
            obj0 = obj;
        end        
        if saveres && mod(iter,5) == 1
            imwrite(kernel./max(kernel(:)),[sdir '/iter' num2str(iter) '_ker.png']);
        end
    end
end
%% 
end








