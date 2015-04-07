function tv = TV(L)
[Lx,Ly] = get_gradient(L);
G = sqrt(Lx.^2+Ly.^2);
tv = sum(G(:));

end

