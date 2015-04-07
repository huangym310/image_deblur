function x = GST(y,lambda,p, J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a generalized threshold shrinkage algorithm for L_p norm minimization
%%  min_{x} = 1/2(x - y)^2 + |x|^p
%% appear in ICCV 2013 Zuo et.al
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Thre = (2*lambda*(1-p))^(1/(2-p)) + lambda*p*(2*lambda*(1-p))^((p-1)/(2-p)); 
n = length(y);
x = zeros(n,1);
for i = 1:n
    if abs(y(i)) < Thre
        x(i) = 0;
    else
        x(i) = abs(y(i));
        for j = 1:J
            x(i) = abs(y(i)) - lambda * p * (x(i)^(p-1));
        end
        x(i) = sign(y(i)) * x(i);
    end
end

