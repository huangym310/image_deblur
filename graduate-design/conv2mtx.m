function [A] = conv2mtx(I,rsize_ker,csize_ker,shape)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% compute matrix A such that conv2(I,kernel) = A * vec(kernel)
if nargin<4
    shape = 'full';
end
[rsize_I csize_I] = size(I);
len_ker = rsize_ker*csize_ker;
    
if strcmpi(shape,'full')
    rsize_I1 = rsize_I + rsize_ker - 1;
    csize_I1 = csize_I + csize_ker - 1;
    A = zeros(rsize_I1*csize_I1,len_ker);

    for i=0:rsize_I1-1
        for j=0:csize_I1-1
            Aij = zeros(1,len_ker);
            for u = max(0,i-rsize_I+1):min(i,rsize_ker-1)
                for v = max(0,j-csize_I+1):min(j,csize_ker-1)
                    ind_ker = v*rsize_ker + u + 1;
                    Aij(ind_ker) = I(i-u+1,j-v+1);
                end
            end
            A(j*rsize_I1+i+1,:) = Aij;
        end
    end
elseif strcmpi(shape,'same')
    A = zeros(rsize_I*csize_I,len_ker);
    rrad_ker = round((rsize_ker-1)/2);
    crad_ker = round((csize_ker-1)/2);
    for i=rrad_ker:rrad_ker+rsize_I-1
        for j=crad_ker:crad_ker+csize_I-1
            Aij = zeros(1,len_ker);
            ind_ij = (j-crad_ker)*rsize_I+i-rrad_ker+1;
            for i1 = max(0,i-rsize_I+1):min(i,rsize_ker-1)
                ind_i = i - i1;
                for j1 = max(0,j-csize_I+1):min(j,csize_ker-1)
                    ind_j = j - j1;
                    ind_ker = j1*rsize_ker + i1 + 1;
                    Aij(ind_ker) = I(ind_i+1,ind_j+1);
                end
                A(ind_ij,:) = Aij;
            end
        end
    end
elseif strcmpi(shape,'valid')
    rsize_I1 = rsize_I - rsize_ker + 1;
    csize_I1 = csize_I - csize_ker + 1;
    A = zeros(rsize_I1*csize_I1,len_ker);
    for i=rsize_ker-1:rsize_I-1
        for j=csize_ker-1:csize_I-1
            Aij = zeros(1,len_ker);
            ind_ij = (j-csize_ker+1)*rsize_I1+i-rsize_ker+1;
            for i1 = 0:rsize_ker-1
                ind_i = i - i1;
                for j1 = 0:csize_ker-1
                    ind_j = j - j1;
                    ind_ker = j1*rsize_ker + i1 + 1;
                    Aij(ind_ker) = I(ind_i+1,ind_j+1);
                end
                A(ind_ij+1,:) = Aij;
            end
        end
    end
else
    error('invalid parameter.');
end

end

