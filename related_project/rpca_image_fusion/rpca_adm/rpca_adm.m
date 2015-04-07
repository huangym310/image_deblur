function [A,E] = rpca_adm(D,lambda)
%   Solves the robust PCA
%      min  ||A||_* + lambda||E||_1
%      s.t.  A+E = D
%   by ADM algorithm.
%   Inputs:
%      D      -- the data matrix, m x n.
%      lambda -- magnitude of 1 norm term
%   Outputs:
%      A -- The estimated of A
%      E -- The estimated of E

%initialize
[m,n] = size(D);
if ~exist('lambda','var')    
    lambda = 1/sqrt(max(m,n));
end
epsilon = 5e-1;
rho = 1.7;
mu_max = 100;
mu = 1.25/svds(D,1);
A = zeros(m,n);
E = zeros(m,n);
Y = zeros(m,n);
iter = 0;
converged = false;
while ~converged
    iter = iter+1; 
    %update E
    E_prev = E;
    E = shrink(D-A+1/mu*Y,lambda/mu);
    %update A  
    A_prev = A;
    [A,rankA] = sv_shrink(D-E+1/mu*Y,1/mu);
    %update Y
    R = D-A-E;
    Y = Y+mu*R; 
    %update mu
    mu = min(mu*rho,mu_max);
    % check convergence
    E_change = E-E_prev;  
    A_change = A-A_prev;
    if max(abs(E_change(:))) < epsilon ...
       && max(abs(R(:))) < epsilon ...
       && max(abs(A_change(:))) < epsilon
        converged = true;
    end
    %display the result   
    disp(['No_SVDs ' num2str(iter) ' r(A) ' num2str(rankA)...
         ' |E|_0 ' num2str(sum(E(:)~=0))...
         ' mu ' num2str(mu)]);     
end% while

end%of function