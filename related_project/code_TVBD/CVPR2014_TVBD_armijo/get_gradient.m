function [Lx,Ly] = get_gradient(L)
% filter: [-1,1]
Lx = L(:,[2:end,end],:)-L;
Ly = L([2:end,end],:,:)-L;

end

