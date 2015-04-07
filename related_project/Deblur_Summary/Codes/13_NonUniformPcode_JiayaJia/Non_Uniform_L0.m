function [retImg, retKernel] = Non_Uniform_L0(filename,ks,wt_l0,wt_deconv,wt_ker)
%Non_Uniform_L0 - Non-Uniform Deblurring via L0 Sparsity
%   [retImg, retKernel] = Non_Uniform_L0(filename,ks,wt_l0,wt_deconv,wt_ker) recovers the latent image "retImg" from
%   the blurred input "filename", with kernel size "ks", image prior weight "wt_l0", and kernel prior weight "wt_ker". 
%   
%   Paras: 
%   @filename  : input image filename; both grayscale and color images are acceptable.
%   @ks    	   : kernel size.
%   @wt_l0     : weight of the L0-norm prior.
%                Typical value in [2e-3, 1e-2];  6e-3 by defalut.                       
%   @wt_deconv : parameter controlling regularization of final non-blind deconvolution. 
%				 The smaller the value is, the smoother the result becomes. 
%                Range [5e2, 8e3]; 2e3 by defalut.   
%   @wt_ker    : weight for the kernel sparsity prior
%            
%   Example
%   ==========
%   [retImg, retKernel] = Non_Uniform_L0('Pantheon.jpg',9,2e-3,4e3,0.1);
%   figure, imshow(retImg), figure, imshow(retKernel,[]);
%
%   ==========
%   The code is created based on the method described in 
%   "Unnatural L0 Sparse Representation for Natural Image Deblurring", Li Xu, Shicheng Zheng, Jiaya Jia, IEEE Conference on 
%   Computer Vision and Pattern Recognition (CVPR) , 2013
%   It is for non-comercial use only.
%  
%   Author: Shicheng Zheng (sczheng@cse.cuhk.edu.hk)
%   Date  : 05/29/2013
%   Version : 1.0 
%   Copyright 2013, The Chinese University of Hong Kong.

end