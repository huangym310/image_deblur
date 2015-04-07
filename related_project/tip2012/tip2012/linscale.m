function R = linscale(I)
%
% R = linscale(I)
% 
% linear rescaling of image I to have values between 0 and 1
%

  R = (I-min(I(:)))/(max(I(:))-min(I(:)));
  