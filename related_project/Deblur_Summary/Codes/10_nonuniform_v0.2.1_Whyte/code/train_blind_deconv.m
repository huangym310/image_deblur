function [x1,x2,error,data_points]=train_blind_deconv(dimensions,ensemble,D,Dp,options)

% Author: James Miskin, adapted by Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

%
% Code from James Miskin, adapted by R.Fergus for use in SIGGRAPH 06
% paper 'Removing Camera Shake from a Single Photograph' by R. Fergus,
% B. Singh, A. Hertzmann, S. T. Roweis and W. T. Freeman
%

% Modified by Oliver Whyte for CVPR 2010 paper:
% "Non-uniform Deblurring for Shaken Images"
% by O. Whyte, J. Sivic, A. Zisserman and J. Ponce

% Original code available at:
% http://www.inference.phy.cam.ac.uk/jwm1003/train_ensemble.tar.gz
% For details of code operation, please see the following publications:
%
% @incollection{Miskin2000,
% author = {James Miskin and David J. C. MacKay},
% title = {Ensemble {L}earning for {B}lind {I}mage {S}eparation and {D}econvolution},
% booktitle = {Adv.~in {I}ndependent {C}omponent {A}nalysis},
% editor = {M. Girolani},
% publisher = {Springer-Verlag},
% year = {2000}
% }
%
% @phdthesis {Miskin2000a,
%        author = "Miskin, James W.",
%        title = {Ensemble Learning for Independent Component Analysis},
%        school = {University of Cambridge},
%        year = {2000},
%        month = {December},
%        url = {http://www.inference.phy.cam.ac.uk/jwm1003/} }
  

%=========================================================================
% This function find the optimal distributions for a deconvolution
% problem in which the data, D, is formed from a convolution of a known
% filter, me, and an unknown source, mx. The data is IxJ pixels, the
% source is MxN pixels and the filter is KxL pixels. The data has a
% binary mask, Dp, which is 1 when the data is observed.
%
% See train_ensemble_demo3.m for example of use.
%=========================================================================

I         = options.I;
J         = options.J;
K         = options.K;
L         = options.L;
M         = options.M;
N         = options.N;     
Kblurry     = options.Kblurry;
Ksharp    = options.Ksharp;
FFT_SHIFT = options.FFT_SHIFT;

nPlanes = dimensions(3,1);
%=========================================================================
% Copy the expectations
%=========================================================================
mmu  = train_ensemble_get(1,dimensions,ensemble.mx);
me   = reshape(train_ensemble_get(2,dimensions,ensemble.mx),K,L);
mx   = reshape(train_ensemble_get(3,dimensions,ensemble.mx),M,N);
mmu2 = train_ensemble_get(1,dimensions,ensemble.mx2);
me2  = reshape(train_ensemble_get(2,dimensions,ensemble.mx2),K,L);
mx2  = reshape(train_ensemble_get(3,dimensions,ensemble.mx2),M,N);

if(FFT_SHIFT)
% compute expanded and shifted convolution kernels
% like fftshift, except we can't use that in the case of a small kernel and
% a large image,
	me = ifftshiftpad(me,[I,J]);
	me2 = ifftshiftpad(me2,[I,J]);
end

%=========================================================================

%=========================================================================
% Evaluate useful FFTs
%=========================================================================
fft_mx =fft2(mx,I,J);
fft_mx2=fft2(mx2,I,J);
fft_mx3=fft2(mx.^2,I,J);
fft_me =fft2(me,I,J);
fft_me2=fft2(me2,I,J);
fft_me3=fft2(me.^2,I,J);

fft_Dp =fft2(Dp,I,J);
mD=real(ifft2(fft_mx.*fft_me));				% mD = mx convolution me
fft_err=fft2(Dp.*(D-mD-mmu),I,J); 			% err = Dp .* (D - mD - mmu)
%=========================================================================

%=========================================================================
%Evaluate the reconstruction error
%=========================================================================
data_points=sum(sum(Dp));
error=sum(sum(Dp.*((D-mD-mmu).^2+real(ifft2(fft_me2.*fft_mx2 - fft_me3.*fft_mx3)))))+data_points*(mmu2-mmu^2);
%=========================================================================

%=========================================================================
%Evaluate the optimal distributions for each parameter
%=========================================================================
mu1=sum(sum(Dp.*(D-mD)));
mu2=data_points;

e1=real(ifft2(fft_err.*conj(fft_mx)));		% e1 = err correlation mx
corr=real(ifft2(fft_Dp.*conj(fft_mx3)));	% corr = Dp correlation mx3

if(FFT_SHIFT)
	e1=e1+me.*corr;	% e1 = e1 + me .* corr
else
	e1(1:K,1:L)=e1(1:K,1:L)+me.*corr(1:K,1:L);
end

e2=real(ifft2(fft_Dp.*conj(fft_mx2)));		% e2 = Dp correlation mx2
x1=real(ifft2(fft_err.*conj(fft_me)));		% x1 = err correlation me
corr=real(ifft2(fft_Dp.*conj(fft_me3)));	% corr = Dp correlation me3

x1(1:M,1:N)=x1(1:M,1:N)+mx.*corr(1:M,1:N);	% x1 = x1 + mx .* corr
x2=real(ifft2(fft_Dp.*conj(fft_me2)));		% x2 = Dp correlation me2

if(FFT_SHIFT)
% undo shift and expansion of convolution kernel
	e1 = fftshiftpad(e1,[K,L]);
	e2 = fftshiftpad(e2,[K,L]);
end

%=========================================================================

%=========================================================================
%Return the optimal distributions
%=========================================================================
x1=[mu1 reshape(e1(1:K,1:L),1,K*L) reshape(x1(1:M,1:N),1,M*N)];
x2=[mu2 reshape(e2(1:K,1:L),1,K*L) reshape(x2(1:M,1:N),1,M*N)];
%=========================================================================


















