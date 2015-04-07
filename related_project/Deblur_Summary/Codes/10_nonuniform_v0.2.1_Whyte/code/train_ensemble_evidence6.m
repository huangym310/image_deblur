function [ensemble,grad] = train_ensemble_evidence6(step_len,dimensions,ensemble,direction,state,D,Dp,options,spatial_blur_mask,not_spatial_image_mask,kernel_mask)
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
% Evaluates evidence for the current ensemble.
%   step_len       - Size of step to take
%   dimensions     - Class description matrix
%   ensemble       - Current ensemble
%   direction      - Direction to search along in ensemble space
%   state          - State variable
%                      1 - Don't update priors or noise
%                      2 - Don't update noise
%                      3 - Update all
%   D-FFT_MODE         - Variables that are passed through to opt_func
%   priors - prior structure with fields pi & gamma
%   spatial_blur_mask - blur mask (vector size of blur) giving weighting on the blur
%   elements
%   not_spatial_image_mask - image mask (vector size of image) giving weighting on variance
%   of image elements - used to alter saturation
%=========================================================================

% Get options
I = options.I;
J = options.J;
K = options.K;
L = options.L;
M = options.M;
N = options.N;
priors = options.priors;
FFT_SHIFT = options.FFT_SHIFT;
Kblurry = options.Kblurry;
theta_grid = options.theta_grid;
non_uniform = options.non_uniform;

Ksharp = Kblurry;
tty=theta_grid{1};
ttx=theta_grid{2};
ttz=theta_grid{3};
if ~exist('kernel_mask')
	kernel_mask = [];
end

numchannels   = dimensions(:,1);
numelements   = dimensions(:,2);
numcomponents = dimensions(:,3);
priortypes    = dimensions(:,4);
priorlock     = dimensions(:,5);
update        = dimensions(:,6);

% ================================================================================
% Update current ensemble by stepping "step_len" along "direction"
% ================================================================================
ensemble.x1     = ensemble.x1         + step_len*direction.x1;
ensemble.x2     = abs(ensemble.x2     + step_len*direction.x2);
ensemble.b_x_2  = abs(ensemble.b_x_2  + step_len*direction.b_x_2);
ensemble.ba_x_2 = abs(ensemble.ba_x_2 + step_len*direction.ba_x_2);
ensemble.pi_x_2 = abs(ensemble.pi_x_2 + step_len*direction.pi_x_2);

% Set to zero any kernel elements contrained to be zero
tmp = train_ensemble_get(2,dimensions,ensemble.x1);
tmp(~kernel_mask) = 1e-16;
ensemble.x1 = train_ensemble_put(2,dimensions,ensemble.x1,tmp);
% Set precision of constrained kernel elements to be large
tmp = train_ensemble_get(2,dimensions,ensemble.x2);
tmp(~kernel_mask) = 1e16*max(tmp(kernel_mask));
ensemble.x2 = train_ensemble_put(2,dimensions,ensemble.x2,tmp);

% ================================================================================
% Initialise new ensemble using direction
% ================================================================================
grad=direction;

if any(isnan(grad.x1)) || any(isnan(grad.x2)), warning('nans in grad'); keyboard; end

% ================================================================================
% Check for valid ensemble, if not, get out
% ================================================================================
if (any(ensemble.x2<0) | any(ensemble.b_x_2<0) | any(ensemble.ba_x_2<0) | any(ensemble.pi_x_2<0))
    ensemble.D_val=Inf;
    grad.x1(:)=NaN;
    grad.x2(:)=NaN;
    grad.b_x_2(:)=NaN;
    grad.ba_x_2(:)=NaN;
    grad.pi_x_2(:)=NaN;
    warning('invalid ensemble')
    return;
end

%Find evidence by summing over the components
% ================================================================================
% For each data type
% ================================================================================
ptr=0;
for c=1:size(dimensions,1)
	% ================================================================================
    % Update the expectations of this data type based on the updated q(x)
	% ================================================================================
    % Get current parameter values of q(x) for this data type
    cx1 =train_ensemble_get(c,dimensions,ensemble.x1);
    cx2 =train_ensemble_get(c,dimensions,ensemble.x2);
    % there are zeros in cx2
    if any(cx2==0)
        disp('zeros in cx2, you''ll get NaNs in opt_c_pi_x_2 and thus NaNs in grad.x1');
    end


    % Evaluate expectations <x> (cmx) and <x^2> (cmx2) under q(x)
    [cHx,cmx,cmx2]=train_ensemble_rectified5(cx1,cx2,priortypes(c));

	% ================================================================================
    % Find optimum values for parameters of prior on each channel of this data type
	% ================================================================================
	% Get current parameter values of prior for this data type
    c_pi_x_2=reshape(ensemble.pi_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
    c_b_x_2=reshape(ensemble.b_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
    c_ba_x_2=reshape(ensemble.ba_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
    c_log_lambda_x=zeros(dimensions(c,1:3));

	% Evaluate log_lambda for each component of this data type's prior
	% I'm not sure what this is, anything to do with the Lagrangian?
    if (numcomponents(c)>1)
		% If mixture has more than one component
        if (priortypes(c)==0 | priortypes(c)==2)
            %Gaussian or Rectified Gaussian prior
            for alpha=1:numcomponents(c)
                for k=1:numchannels(c)
                    c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/ c_pi_x_2(k,alpha)+0.5*log(c_ba_x_2(k,alpha))-0.25/c_b_x_2(k,alpha)-0.5*cmx2(k,:)*c_ba_x_2(k,alpha));
                end
            end
        elseif (priortypes(c)==1)
            %Exponential prior
            for alpha=1:numcomponents(c)
                if isempty(spatial_blur_mask)
					disp('isempty(spatial_blur_mask)')
                    keyboard
                    for k=1:numchannels(c)
                        c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+log(c_ba_x_2(k,alpha))-0.5/c_b_x_2(k,alpha)-cmx(k,:)*c_ba_x_2(k,alpha)) + spatial_blur_mask(alpha,:);
                    end
                else
                    for k=1:numchannels(c)
                        c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+log(c_ba_x_2(k,alpha))-0.5/c_b_x_2(k,alpha)-cmx(k,:)*c_ba_x_2(k,alpha));
                    end
                end
            end
            if nnz(isnan(c_log_lambda_x))
                disp('c_log_lambda_x has nans in it');
            end
        elseif (priortypes(c)==3)
            %Discrete prior
        elseif  (priortypes(c)==4)
            %Exponential prior
            for alpha=1:numcomponents(c)
                for k=1:numchannels(c)
                    c_log_lambda_x(k,:,alpha)=log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+0.5*log(c_ba_x_2(k,alpha))-0.25/c_b_x_2(k,alpha)-0.5*abs(cmx(k,:))*c_ba_x_2(k,alpha);
                end
            end
        end

        % Now normalise c_log_lambda_x
		% Take max over all components of the prior
        max_c_log_lambda_x=max(c_log_lambda_x,[],3);
        for alpha=1:numcomponents(c)
            c_log_lambda_x(:,:,alpha)=c_log_lambda_x(:,:,alpha)-max_c_log_lambda_x;
        end
		% Sum lambdas over all components of the prior and re-take logs
        log_sum_lambda_x=log(sum(exp(c_log_lambda_x),3));
        for alpha=1:numcomponents(c)
            c_log_lambda_x(:,:,alpha)=c_log_lambda_x(:,:,alpha)-log_sum_lambda_x;
        end
    end

    % Using log_lambda, compute optimum parameter values for priors for each channel of this data type
    if (priortypes(c)==0 | priortypes(c)==2)
        %(Rectified) Gaussian prior for P(x)
        opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c))/2;
        opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c));
        opt_c_ba_x_2=zeros(numchannels(c),numcomponents(c));
        for alpha=1:numcomponents(c)
            opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx2,2)/2);
        end
    elseif (priortypes(c)==1)
        %Exponential prior for P(x)
        opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c));
        opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c));
        opt_c_ba_x_2=zeros(numchannels(c),numcomponents(c));
        for alpha=1:numcomponents(c)
            opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx,2));
        end
    elseif ((priortypes(c)==3) | (priortypes(c)==4))
        %Discrete prior (not learning any parameters of prior)
        opt_c_b_x_2=ensemble.b_x*ones(numchannels(c),numcomponents(c));
        opt_c_ba_x_2=ensemble.b_x/ensemble.a_x*ones(numchannels(c),numcomponents(c));
        opt_c_pi_x_2=ensemble.pi_x*ones(numchannels(c),numcomponents(c));
    elseif (priortypes(c)==5)
        %Laplacian prior for P(x)
        opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c));
        opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),numchannels(c),numcomponents(c));
        opt_c_ba_x_2=zeros(numchannels(c),numcomponents(c));
        for alpha=1:numcomponents(c)
            opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx,2));
        end
    end

	% ================================================================================
    % Manual over-ride for priors
	% ================================================================================
    if (priorlock(c)>0)
		% if priorlock is set for this data type
        if (c==3) %% Image prior
            opt_c_pi_x_2 = priors.pi(1:numchannels(c),:) * 1e3;
            opt_c_ba_x_2 = priors.gamma(1:numchannels(c),:);
            opt_c_b_x_2 = ones(numchannels(c),numcomponents(c)) * 1e-3;
        elseif (c==2) %% Blur prior
            opt_c_ba_x_2 = [5.1143e3 5.0064e+03 173.8885 50.6538 ];
            opt_c_b_x_2  = [787.8988 201.7349   236.1948 143.1756];
        else
            error foo
        end
    else
        % Leave alone
    end
    if any(isnan(opt_c_pi_x_2))
        disp('nans in opt_c_pi_x_2');
    end

	% ================================================================================
    % Having found optimum parameter values for priors, begin calculating new parameters for q(x) under these priors, by adding the data-independent, parameter-derived terms
	% ================================================================================
    %Optimum Q(x)
    opt_cx1=zeros(size(cx1));
    opt_cx2=zeros(size(cx2));
    if (priortypes(c)==0 | priortypes(c)==2)
		% if prior is Gaussian or rectified Gaussian
        for alpha=1:numcomponents(c) % for each mixture componenent
            for k=1:numchannels(c) % for each hidden image
                opt_cx2(k,:)=opt_cx2(k,:)+c_ba_x_2(k,alpha)*exp(c_log_lambda_x(k,:,alpha));
            end
        end
    elseif ((priortypes(c)==1)| priortypes(c)==4)
		% if prior is exponential or discrete
        for alpha=1:numcomponents(c)
            for k=1:numchannels(c)
                opt_cx1(k,:)=opt_cx1(k,:)-c_ba_x_2(k,alpha)*exp(c_log_lambda_x(k,:,alpha));
            end
        end
    end

	% ================================================================================
    % Evaluate the KL divergence for this data type under current ensemble
	% ================================================================================
    ensemble.D_x(sum(dimensions(1:c-1,1))+1:sum(dimensions(1:c,1)))=sum(cHx,2)+...
        sum(gammaln(ensemble.b_x)-ensemble.b_x*log(ensemble.a_x)-gammaln(c_b_x_2)+c_b_x_2.*log(c_b_x_2./c_ba_x_2)+...
        (c_b_x_2-opt_c_b_x_2).*(log(c_ba_x_2)-0.5./c_b_x_2)+(opt_c_b_x_2./opt_c_ba_x_2-c_b_x_2./c_ba_x_2).*c_ba_x_2,2)+...
        sum(gammaln(ensemble.pi_x)-gammaln(c_pi_x_2)+(c_pi_x_2-opt_c_pi_x_2).*(log(c_pi_x_2)-0.5./c_pi_x_2),2)+...
        (-gammaln(numcomponents(c)*ensemble.pi_x)+gammaln(sum(c_pi_x_2,2))+sum(c_pi_x_2-opt_c_pi_x_2,2).*(-log(sum(c_pi_x_2,2))+0.5./sum(c_pi_x_2,2)))+...
        sum(sum(c_log_lambda_x.*exp(c_log_lambda_x),2),3);

	% ================================================================================
    % Store updated prior parameters and expectations under current distribution q(x)
	% ================================================================================
    ensemble.log_lambda_x=train_ensemble_put_lambda(c,dimensions,ensemble.log_lambda_x,c_log_lambda_x);
    ensemble.mx=train_ensemble_put(c,dimensions,ensemble.mx,cmx);
    ensemble.mx2=train_ensemble_put(c,dimensions,ensemble.mx2,cmx2);
    ensemble.opt_ba_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(opt_c_ba_x_2,1,numchannels(c)*numcomponents(c));
    ensemble.opt_b_x_2( ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(opt_c_b_x_2 ,1,numchannels(c)*numcomponents(c));
    ensemble.opt_pi_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(opt_c_pi_x_2,1,numchannels(c)*numcomponents(c));

	% ================================================================================
	%  Put data-independent part of new parameters of q into new ensemble
	% ================================================================================
    grad.x1=train_ensemble_put(c,dimensions,grad.x1,opt_cx1);
    grad.x2=train_ensemble_put(c,dimensions,grad.x2,opt_cx2);

	% ================================================================================
    % Increment the pointer to point to start of priors for next data-type
	% ================================================================================
    ptr=ptr+numchannels(c)*numcomponents(c);

	if any(isnan(grad.x1)) || any(isnan(grad.x2))
		error('NaNs in grad.x1 or grad.x2 in train_ensemble_evidence6');
	end
end

% ================================================================================
% Only set the gradient for the priors if in state 2 or 3
% ================================================================================
if (state>=2)
	% set prior parameter values of new ensemble to optimal values for current ensemble
    grad.pi_x_2=ensemble.opt_pi_x_2;
    grad.b_x_2 =ensemble.opt_b_x_2;
    grad.ba_x_2=ensemble.opt_ba_x_2;
else
	% otherwise, set prior parameter values of new ensemble to those of current ensemble
    grad.pi_x_2=ensemble.pi_x_2;
    grad.b_x_2 =ensemble.b_x_2;
    grad.ba_x_2=ensemble.ba_x_2;
end

%Q(gamma)
% ================================================================================
% Finish calculating the parameters for the new q(x), adding the data-dependent terms calculated in train_blind_deconv.
% ================================================================================
if non_uniform
	% which train_blind_deconv2
    [dx1,dx2,rerror,data_points]=train_blind_deconv2(dimensions,ensemble,D,Dp,options,ttx,tty,ttz);
else
	[dx1,dx2,rerror,data_points]=train_blind_deconv(dimensions,ensemble,D,Dp,options);
end

% ================================================================================
% Update parameters of q(beta_sigma) (noise parameter) for new ensemble
% ================================================================================
ensemble.b_sigma_2=ensemble.b_sigma+data_points/2; % Same for current and new ensemble
ensemble.opt_ba_sigma_2=ensemble.b_sigma_2/(ensemble.a_sigma+rerror/2);
% Set current ensemble's noise parameters to the optimum if in the final state
if (state==3)
    ensemble.ba_sigma_2=ensemble.opt_ba_sigma_2;
end

% ================================================================================
% Add data-dependent terms of parameters of q(x) to new ensemble
% ================================================================================
grad.x1 = grad.x1+ensemble.ba_sigma_2*dx1;
grad.x2 = grad.x2+ensemble.ba_sigma_2*dx2;
if any(isnan(grad.x1)) || any(isnan(grad.x2)), warning('nans in grad'); keyboard; end
% ================================================================================
% Evaluate the KL divergence of q(beta_sigma) under current ensemble
% ================================================================================
ensemble.D_x(sum(numchannels)+1)=gammaln(ensemble.b_sigma)-...
    gammaln(ensemble.b_sigma_2)-...
    ensemble.b_sigma*log(ensemble.a_sigma)+...
    ensemble.b_sigma_2*log(ensemble.b_sigma_2/ensemble.ba_sigma_2)+...
    (ensemble.b_sigma_2/ensemble.opt_ba_sigma_2-ensemble.b_sigma_2/ensemble.ba_sigma_2)*ensemble.ba_sigma_2+...
    data_points*log(2*pi)/2;

% ================================================================================
% Normalise KL divergences from log_e to bits per data point
% ================================================================================
ensemble.D_x=ensemble.D_x/data_points*log2(exp(1));
if (isnan(ensemble.D_val))
    ensemble.D_val = inf;
end

% Evaluate the total KL divergence
ensemble.D_val = sum(ensemble.D_x);

% ================================================================================
% Subtract the current ensemble from the new one to find the 'gradient'
% ================================================================================
grad.x1     = grad.x1-ensemble.x1;
grad.x2     = grad.x2-ensemble.x2;
grad.b_x_2  = grad.b_x_2-ensemble.b_x_2;
grad.ba_x_2 = grad.ba_x_2-ensemble.ba_x_2;
grad.pi_x_2 = grad.pi_x_2-ensemble.pi_x_2;

% ================================================================================
% Set grad to zero for any kernel elements that are clamped
% ================================================================================
tmp               = train_ensemble_get(2,dimensions,grad.x1);
tmp(~kernel_mask) = 0;
grad.x1           = train_ensemble_put(2,dimensions,grad.x1,tmp);
tmp               = train_ensemble_get(2,dimensions,grad.x2);
tmp(~kernel_mask) = 0;
grad.x2           = train_ensemble_put(2,dimensions,grad.x2,tmp);

% ================================================================================
% Set grad to zero for any classes that are clamped; don't train classes if options specify clamping
% ================================================================================
ptr=0;
for c=1:size(dimensions,1)
    if (~update(c))
        %%% Comment out line below & set dimensions(c,6) to 0 for MAP
        %%% estimate
        grad.x1(ptr+1:ptr+numchannels(c)*numelements(c)) = 0;
        grad.x2(ptr+1:ptr+numchannels(c)*numelements(c)) = 0;
    end
    ptr = ptr+numchannels(c)*numelements(c);
end
