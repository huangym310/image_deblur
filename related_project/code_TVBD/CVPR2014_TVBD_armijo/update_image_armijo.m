function [L,learning_rate] = update_image_armijo(B,K,L,eta,lambda,learning_rate)
% update image using gradient descent with Armijo rule
sigma = 0.1;
beta = 0.5;

obj = objective_fun(B,K,L,eta,lambda);
grad = gradient_fun(B,K,L,eta,lambda);
neg_grad = -grad;     
obj_new = objective_fun(B,K,L+learning_rate*neg_grad,eta,lambda);
iter = 0;
while (obj_new-obj)/learning_rate > sigma*grad'*neg_grad
    learning_rate = learning_rate*beta;
    obj_new = objective_fun(B,K,L+learning_rate*neg_grad,eta,lambda);
    iter = iter+1;
end%while
L = L+learning_rate*neg_grad;
% fprintf('Update image for %d times!\n',iter);
   
end%function

function obj = objective_fun(B,K,L,eta,lambda)
n = size(B,3);
obj = lambda*TV(L);
for ii = 1:n
    E = B(:,:,ii)-conv2(L,K(:,:,ii),'valid');
    obj = obj+eta(ii)/2*E(:)'*E(:);    
end%ii

end%function

function dL = gradient_fun(B,K,L,eta,lambda)
n = size(B,3);
dL = lambda*gradTVcc(L);
for ii = 1:n
    dL = dL+eta(ii)*conv2(conv2(L,K(:,:,ii),'valid')-B(:,:,ii),rot90(K(:,:,ii),2),'full');
end

end
