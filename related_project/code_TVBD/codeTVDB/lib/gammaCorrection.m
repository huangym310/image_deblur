function [Ic] = gammaCorrection(I,gamma)
    Ic = I.^(gamma);
end