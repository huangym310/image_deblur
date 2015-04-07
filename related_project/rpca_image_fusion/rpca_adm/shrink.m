function SW = shrink(W,epsilon)
% soft-thresholding (shrinkage) operator
% solve the following optimization problem:
% SW = arg_X  min epsilon||X||_1 + 1/2||X-W||_F^2

SW = sign(W).*max(abs(W)-epsilon,0);

end%of function