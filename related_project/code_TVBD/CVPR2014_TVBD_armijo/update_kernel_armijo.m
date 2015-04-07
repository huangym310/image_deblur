function [K,eta_K] = update_kernel_armijo(B,K,L,gamma,eta_K)
% compute kernel prior
[Lx,Ly] = get_gradient(L);
filt = ones(5)/5^2;
Lx = conv2(Lx,filt,'same');
Ly = conv2(Ly,filt,'same');
G = sqrt(Lx.^2+Ly.^2);
PK = (K>0);
% update kernel using gradient descent with Armijo rule
n = size(B,3);
for ii = 1:n 
    [K(:,:,ii),eta_K(ii)] = update_kernel_single(B(:,:,ii),K(:,:,ii),L,...
                                                 G,PK(:,:,ii),gamma,eta_K(ii));
end

end%function

function [K,eta_K] = update_kernel_single(B,K,L,G,PKi,gamma,eta_K)
sigma = 0.1;
beta = 0.5;
obj = objective_fun(B,K,L,G,PKi,gamma);
grad = gradient_fun(B,K,L,G,PKi,gamma);
neg_grad = -grad;     
obj_new = objective_fun(B,K+eta_K*neg_grad,L,G,PKi,gamma);
iter = 0;
while (obj_new-obj)/eta_K > sigma*grad'*neg_grad
    eta_K = eta_K*beta;
    obj_new = objective_fun(B,K+eta_K*neg_grad,L,G,PKi,gamma);
    iter = iter +1;
end%while
K = K+eta_K*neg_grad;
% fprintf('Update kernel for %d times!\n',iter);
   
end%function

function obj = objective_fun(B,K,L,G,PKi,gamma)
E = B-conv2(L,K,'valid');
T = conv2(G,K.*PKi,'valid');
obj = 1/2*E(:)'*E(:)+gamma/2*T(:)'*T(:); 
end%function

function dK = gradient_fun(B,K,L,G,PKi,gamma)
dK = conv2(rot90(L,2),conv2(L,K,'valid')-B,'valid')...
    +gamma*PKi.*conv2(rot90(G,2),conv2(G,K,'valid'),'valid');
end