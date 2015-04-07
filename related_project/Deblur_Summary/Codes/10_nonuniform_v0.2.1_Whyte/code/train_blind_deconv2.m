function [x1,x2,rerror,data_points]=train_blind_deconv2(dimensions,ensemble,D,Dp,options,ttx,tty,ttz)
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

I           = options.I;
J           = options.J;
K           = options.K;
L           = options.L;
M           = options.M;
N           = options.N;     
Kblurry     = options.Kblurry;
Ksharp      = options.Ksharp;
NUM_THREADS = options.NUM_THREADS;

%=========================================================================
% Copy the expectations
%=========================================================================
mmu  = train_ensemble_get(1,dimensions,ensemble.mx);
me   = reshape(train_ensemble_get(2,dimensions,ensemble.mx),K,L);
mx   = reshape(train_ensemble_get(3,dimensions,ensemble.mx),M,N);
mmu2 = train_ensemble_get(1,dimensions,ensemble.mx2);
me2  = reshape(train_ensemble_get(2,dimensions,ensemble.mx2),K,L);
mx2  = reshape(train_ensemble_get(3,dimensions,ensemble.mx2),M,N);
%=========================================================================

theta_grid = [ttx(:),tty(:),ttz(:)]';
clamp_edges_to_zero = 0;
non_uniform = 1;

% Test to check both versions give same outputs
% [x1a,x2a,e1a,e2a,mu1a,mu2a,errora] = train_blind_deconv_mex(mx,mx2,me,me2,D,Ksharp,Kblurry,-theta_grid,Dp,mmu,mmu2);
% [x1b,x2b,e1b,e2b,mu1b,mu2b,errorb] = train_blind_deconv_mex_threads(mx,mx2,me,me2,D,Ksharp,Kblurry,-theta_grid,Dp,mmu,mmu2,NUM_THREADS);
% fprintf('x1: %f, x2: %f, e1: %f, e2: %f, mu1: %f, mu2: %f, error: %f\n', ...
% 		sqrt(sum((x1a(:)-x1b(:)).^2)), sqrt(sum((x2a(:)-x2b(:)).^2)), ...
% 		sqrt(sum((e1a(:)-e1b(:)).^2)), sqrt(sum((e2a(:)-e2b(:)).^2)), ...
% 		sqrt(sum((mu1a(:)-mu1b(:)).^2)), sqrt(sum((mu2a(:)-mu2b(:)).^2)), sqrt(sum((errora(:)-errorb(:)).^2)));

if NUM_THREADS > 1
	% Requires pthreads
	[x1,x2,e1,e2,mu1,mu2,rerror] = train_blind_deconv_mex_threads(mx,mx2,me,me2,D,Ksharp,Kblurry,-theta_grid,double(Dp),mmu,mmu2,NUM_THREADS);
else
	% No pthreads
	[x1,x2,e1,e2,mu1,mu2,rerror] = train_blind_deconv_mex(mx,mx2,me,me2,D,Ksharp,Kblurry,-theta_grid,Dp,mmu,mmu2);
end
data_points = mu2;

%=========================================================================
%Return the optimal distributions
%=========================================================================
x1=[mu1 reshape(e1(1:K,1:L),1,K*L) reshape(x1(1:M,1:N),1,M*N)];
x2=[mu2 reshape(e2(1:K,1:L),1,K*L) reshape(x2(1:M,1:N),1,M*N)];
%=========================================================================

