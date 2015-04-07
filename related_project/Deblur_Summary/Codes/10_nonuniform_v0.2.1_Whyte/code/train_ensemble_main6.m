function [ensemble,D_log,gamma_log]=train_ensemble_main6(dimensions,initial_x1,initial_x2,options,D,Dp,spatial_blur_mask,not_spatial_image_mask,kernel_mask)
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
% Performs ensemble learning to model a data set.
%   dimensions     - Class description matrix
%   initial_x1,x2  - Initial parameters of the ensemble, can be []
%   text           - Text that is displayed at each iteration
%   options        - Options to control the training
%   D-FFT_MODE         - Variables that are passed through to opt_func
%=========================================================================
% Performs ensemble fitting of posterior over parameters by a separable
% mixture distribution. The problem is described by the dimensions array
% which stores a row for each type of data in the problem. Each row has
% the following columns (for row i)
%  1 - Number of rows of x_i (each row has a different prior)
%  2 - Number of columns of x_i
%  3 - Number of components to use for the mixture prior
%  4 - Type of prior to use for x_i
%        0 - Gaussian
%        1 - Exponential
%        2 - Rectified Gaussian
%
% When training, opt_func is called to evaluate the optimal
% distributions. It is called as
%   opt_func(dimensions,ensemble,D,Dp,I,...)
% The routine must return the number of data points, the error in the
% reconstruction and the optimum values for the parameters of the
% distributions (without including the prior). See the example models to
% see how this works. The priors that have been implemented so far result
% in a posterior is of the form
%    log Q(x) = -x2*x/2+x1*x + constant
% so the model file returns the optimum x1,x2 parameters given the
% data. The library routines add in the terms from the priors.
%
% options controls the training
%   options(1) - Convergence criteria, the algorithm is assumed to be
%                converging if the improvement at each iteration is less
%                than the criteria.
%   options(2) - non-zero if the cost function along the search direction
%                is to be displayed at each iteration (slow).
%   options(3) - Initial noise variance to be assumed.
%   options(4) - Non-zero if priors are to be reinitialised whenever the
%                noise changes
%   options(5) - Non-zero if any components that are switched off (ie
%                tend to prior) are to be re-randomised when noise
%                changes, this may be an advantage if the noise is
%                initially assumed to be large.
%=========================================================================

% reconstruct observed gradients into intensity image
% Drecon = reconsEdge3(D);

%Create array to store training logs
Niter     = options.num_iters;
D_log     = NaN*ones(2,Niter);
gamma_log = NaN*ones(1,Niter);

if ~exist('spatial_blur_mask')
    spatial_blur_mask = [];
end

if ~exist('not_spatial_image_mask')
    not_spatial_image_mask = [];
end

if ~exist('kernel_mask')
	kernel_mask = [];
end

I                        = options.I;
J                        = options.J;
K                        = options.K;
L                        = options.L;
M                        = options.M;
N                        = options.N;
priors                   = options.priors;
text                     = options.text;
FFT_SHIFT                = options.FFT_SHIFT;
Kblurry                  = options.Kblurry;
theta_grid               = options.theta_grid;
non_uniform              = options.non_uniform;
SAVE_EACH_ITERATION = options.SAVE_EACH_ITERATION;
if isfield(options,'OUTPUT_DIRECTORY'), OUTPUT_DIRECTORY = options.OUTPUT_DIRECTORY; else OUTPUT_DIRECTORY = []; end
if isfield(options,'OUTPUT_FILENAME'), OUTPUT_FILENAME = options.OUTPUT_FILENAME; else OUTPUT_FILENAME = []; end
if isfield(options,'s'), s = options.s; else s = []; end
	
% delete any per-iteration images from previous runs
delete(fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,sprintf('s%02d_it*.jpg',s)));

% ================================================================================
% Copy dimensions
% ================================================================================
numchannels   = dimensions(:,1);
numelements   = dimensions(:,2);
numcomponents = dimensions(:,3);
priortypes    = dimensions(:,4);
priorlock     = dimensions(:,5);

% ================================================================================
% Create the ensemble parameters
% ================================================================================
% a / b are 1st and 2nd hyperparameters for distributions, eg mean / variance for a Gaussian.
ensemble=struct('x1' ,1e4*randn(1,numchannels'*numelements).*ceil(rand(1,numchannels'*numelements)*2),...
    'x2',1e4*ones(1,numchannels'*numelements),...
    'mx' ,zeros(1,numchannels'*numelements),...
    'mx2',zeros(1,numchannels'*numelements),...
    'log_lambda_x',zeros(1,sum(numchannels.*numelements.*numcomponents)),...
    'pi_x',1,...
    'pi_x_2',ones(1,numchannels'*numcomponents),...
    'opt_pi_x_2',ones(1,numchannels'*numcomponents),...
    'a_x',1e-3,...
    'ba_x_2' ,ones(1,numchannels'*numcomponents),...
    'opt_ba_x_2' ,ones(1,numchannels'*numcomponents),...
    'b_x',1e-3,...
    'b_x_2' ,ones(1,numchannels'*numcomponents),...
    'opt_b_x_2' ,ones(1,numchannels'*numcomponents),...
    'a_sigma',1e-3,...
    'ba_sigma_2',0,...
    'opt_ba_sigma_2',0,...
    'b_sigma',1e-3,...
    'b_sigma_2',0,...
    'D_val',0,...
    'D_x',zeros(sum(numchannels)+1,1));

direction=struct('x1' ,zeros(1,numchannels'*numelements),...
    'x2' ,zeros(1,numchannels'*numelements),...
    'pi_x_2',zeros(1,numchannels'*numcomponents),...
    'ba_x_2' ,zeros(1,numchannels'*numcomponents),...
    'b_x_2' ,zeros(1,numchannels'*numcomponents));

% ================================================================================
% Initialise ensemble if specified
% ================================================================================
if (~isempty(initial_x1))
    ensemble.x1    =initial_x1;
end
if (~isempty(initial_x2))
    ensemble.x2    =initial_x2;
end



% ================================================================================
% Make sure no parameters are set to zero initially
% ================================================================================
ensemble.x1 = ensemble.x1+1e-16*(ensemble.x1==0);
ensemble.x2 = ensemble.x2+1e-16*(ensemble.x2==0);

% Details of the state
% State alternates between 1 and 2, switching every time the KL divergence converges but the noise does not. When both have stabilised, then state is set to 3 for the final iteration
% state=1:	Components of the priors may be reset if they have coalesced or converged to zero.
% 			Channels of the ensemble may be randomised if they have converged to the prior (when q, our approximation to the posterior, converges to the prior, it's probably not getting any information from the data, so we may as well kick start it to have another go at extracting some information from the observations). 
% 			train_ensemble_evidence6 will only optimise the latent variables in the line search, not the priors or the noise parameters
% 			if the KL divergence has converged in a state=1 iteration, the noise will be set directly to the optimum for the final ensemble
% state=2:	
% 			train_ensemble_evidence6 will optimise the latent variables and the priors in the line search
% state=3:	final iteration before exit. Do not reset / randomise any prior components or expectations.
% 			train_ensemble_evidence6 will optimise the latent variables, the priors and the noise parameter in the line search
state            = 1;
last_change_iter = 0;
alpha            = 1.0;
beta             = 0.9;

% ================================================================================
% Copy options
% ================================================================================
converge_criteria      = options.converge_criteria;
plot_linesearch        = options.plot_linesearch;
DISPLAY_EACH_ITERATION = options.DISPLAY_EACH_ITERATION && ~isdeployed;
ensemble.ba_sigma_2    = options.init_noise_var.^-2;
restart_priors         = options.restart_priors;
restart_switched_off   = options.restart_switched_off;

% ================================================================================
% Iterate algorithm to improve the approximating ensemble
% ================================================================================
oD_val=NaN;
for iter=1:Niter
	% ================================================================================
    % Re-evaluate after a state change. Reset step_len to 1.0
	% ================================================================================
    if (iter==last_change_iter+1)
        % if state < 3, we are not in the last iteration.
        % evaluate evidence of current ensemble, adjust priors to their
        % optimal parameters, and reset any components of priors that have 
        % converged to the same point. Get search direction for whole 
        % ensemble.
        if (state<3)
			% If state < 3, we're not in the last iteration
			% ================================================================================
			% Update parameter values of priors using current ensemble
			% ================================================================================
            % Evaluate the current model before updating the priors, hence step_len == 0
            step_len=0.0;
            % Find direction towards optimal ensemble from current one, and optimal prior parameter values
            [ensemble,grad] = train_ensemble_evidence6(step_len,dimensions,ensemble,direction,state,D,Dp,options,spatial_blur_mask,not_spatial_image_mask,kernel_mask);
            direction=grad;

            % Set the priors to their optimal values -- should do nothing, due to priorlock
            ensemble.pi_x_2 = ensemble.opt_pi_x_2;
            ensemble.b_x_2  = ensemble.opt_b_x_2;
            ensemble.ba_x_2 = ensemble.opt_ba_x_2;

            % Re-initialise components of priors where a component has zero weight or two components have tended to the same point
            ptr=0;
            for c=1:size(dimensions,1)
                if (numcomponents(c)>1)
	                % Only need to do this if the prior that have more than one component
                    c_pi_x_2 = reshape(ensemble.pi_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
                    c_b_x_2  = reshape(ensemble.b_x_2( ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
                    c_ba_x_2 = reshape(ensemble.ba_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)),numchannels(c),numcomponents(c));
                    % for each channel
                    for k=1:numchannels(c)
                        %Check for each distribution type for coallescence
                        if (priortypes(c)==0 | priortypes(c)==1 | priortypes(c)==2)
                            %Gaussian or Exponential
                            sorted_scales=sort(c_ba_x_2(k,:));
                            if (restart_priors | any(c_pi_x_2(k,:)<ensemble.pi_x+1/numcomponents(c)) | any(sorted_scales(2:end)<1.5*sorted_scales(1:end-1)))
                                %If any component has approximately zero weight or if any
                                %two have approximately the same scale, they should be
                                %restarted
                                mean_scale=sum(c_b_x_2(k,:)./c_ba_x_2(k,:))/sum(c_b_x_2(k,:));
                                c_pi_x_2(k,:) = ensemble.pi_x+numelements(c)/numcomponents(c);
                                c_b_x_2(k,:)  = ensemble.b_x +numelements(c)/numcomponents(c);
                                c_ba_x_2(k,:)=c_b_x_2(k,:)./(ensemble.a_x+0.5*(1:numcomponents(c))*mean_scale*numelements(c)/numcomponents(c));
                            end
                        elseif (priortypes(c)==3 | priortypes(c)==4)
                            %Discrete so no prior properties
                        end
                    end

                    %Store the new, restarted prior parameters
                    ensemble.pi_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(c_pi_x_2,1,numchannels(c)*numcomponents(c));
                    ensemble.b_x_2( ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(c_b_x_2 ,1,numchannels(c)*numcomponents(c));
                    ensemble.ba_x_2(ptr+1:ptr+numchannels(c)*numcomponents(c)) = reshape(c_ba_x_2,1,numchannels(c)*numcomponents(c));
                end

                % Increment the pointer by the number of channels in
                % this data type times the number of components in the data
                % type's prior, since each channel in the data type
                % gets its own parameters for the prior, even though the
                % form of the prior is shared
                ptr=ptr+numchannels(c)*numcomponents(c);
            end
        end

		% ================================================================================
        % Re-evaluate evidence for current ensemble and get search direction
		% ================================================================================
		step_len=0.0;
        [ensemble,grad]=train_ensemble_evidence6(step_len,dimensions,ensemble,direction,state,D,Dp,options,spatial_blur_mask,not_spatial_image_mask,kernel_mask);

		% ================================================================================
        % Set search direction to the direction of the optimal ensemble under the current distribution, which is given by grad
		% ================================================================================
        direction.x1=alpha*grad.x1;
        direction.x2=alpha*grad.x2;
        direction.b_x_2=alpha*grad.b_x_2;
        direction.ba_x_2=alpha*grad.ba_x_2;
        direction.pi_x_2=alpha*grad.pi_x_2;
        step_len=1.0;
    end

	% ================================================================================
    % Evaluate evidence for intermediate ensembles at various step lengths along search direction, plot KL divergences
	% ================================================================================
    if (plot_linesearch)
        ostep_len=step_len;
        x_vals=(-2:0.2:2)*step_len;
        f_vals=zeros(size(x_vals));
        for i=1:prod(size(x_vals))
            step_len=x_vals(i);
            [pensemble,pgrad]=train_ensemble_evidence6(step_len,dimensions,ensemble,direction,state,D,Dp,options,spatial_blur_mask,not_spatial_image_mask,kernel_mask);
            f_vals(i)=pensemble.D_val; % Store total KL divergence of pensemble
        end
        figure(10); clf;
        plot(x_vals,f_vals)
		title(['Iter ',num2str(iter)])
        % v=axis;
        % line([ostep_len ostep_len],v(3:4))
		drawnow
        step_len=ostep_len;
    end

	% ================================================================================
    % Find minimum along search direction for a step_len between 0 and 2
	% ================================================================================
	% Evaluate evidence for tensemble (=ensemble+direction) and test for improvement
	[tensemble,tgrad] = train_ensemble_evidence6(step_len, dimensions, ensemble, direction, state, D, Dp, options, spatial_blur_mask, not_spatial_image_mask, kernel_mask);
	success=(tensemble.D_val<ensemble.D_val);
	tmpcount = 0;
	% If no improvement, do a kind of line search along the search direction
	if (tensemble.D_val>ensemble.D_val)
		direction.x1=alpha*grad.x1;
		direction.x2=alpha*grad.x2;
		direction.b_x_2=alpha*grad.b_x_2;
		direction.ba_x_2=alpha*grad.ba_x_2;
		direction.pi_x_2=alpha*grad.pi_x_2;
		step_len=2*step_len;
		while (tensemble.D_val>ensemble.D_val+converge_criteria/1e4)
			tmpcount = tmpcount + 1;
			step_len=0.5*step_len;
			[tensemble,tgrad] = train_ensemble_evidence6(step_len, dimensions, ensemble, direction, state, D, Dp, options, spatial_blur_mask, not_spatial_image_mask, kernel_mask);
		end
		if DISPLAY_EACH_ITERATION && iter>1,figure(1);subplot(1,3,3);plot(iter,tensemble.D_val,'ro');end
	else
		if DISPLAY_EACH_ITERATION && iter>1,figure(1);subplot(1,3,3);plot(iter,tensemble.D_val,'go');end
	end

	% ================================================================================
	% ?
	% ================================================================================
	direction.x1     = alpha*tgrad.x1+beta*direction.x1;
	direction.x2     = alpha*tgrad.x2+beta*direction.x2;
	direction.b_x_2  = alpha*tgrad.b_x_2+beta*direction.b_x_2;
	direction.ba_x_2 = alpha*tgrad.ba_x_2+beta*direction.ba_x_2;
	direction.pi_x_2 = alpha*tgrad.pi_x_2+beta*direction.pi_x_2;
	step_len         = min(1,1.1*step_len);

	% Find change in KL divergence from previous iteration
    dD_val = tensemble.D_val-oD_val;

    % Accept the new ensemble and replace current ensemble
    ensemble = tensemble;
    grad     = tgrad;

    % Record KL divergence and noise parameter into log
    D_log(1,iter)=ensemble.D_val;
    gamma_log(1,iter)=ensemble.ba_sigma_2;

	% Display current state
	imtmp = reconsEdge3(reshape(train_ensemble_get(3,dimensions,ensemble.mx),M,N));
	if DISPLAY_EACH_ITERATION
		figure(1)
		% Show image
		subplot(1,3,1);
		imagesc(imtmp);
		title(['Iter ',num2str(iter)])
		axis image; axis off; colormap gray;
		% Show kernel
		subplot(1,3,2);
		if non_uniform
			plot_nonuni_kernel(train_ensemble_get(2,dimensions,ensemble.mx),theta_grid,1,0,0,0,1);
		else
			imagesc(reshape(train_ensemble_get(2,dimensions,ensemble.mx),K,L)); colormap gray; axis off; axis image;
		end
		% Plot cost function
		subplot(1,3,3)
		if iter > 1
			hold on; 
			plot([iter-1:iter],D_log(1,iter-1:iter));
			title(['Iter ',num2str(iter)])
			if iter==last_change_iter+1
				plot(last_change_iter,D_log(1,last_change_iter),'rx');
			end
		else
			cla; hold on; 
		end
		drawnow
	end
    
	if SAVE_EACH_ITERATION && ~isempty(OUTPUT_DIRECTORY) && ~isempty(OUTPUT_FILENAME)
		imwrite((imtmp-min(imtmp(:)))/(max(imtmp(:))-min(imtmp(:))),fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,sprintf('s%02d_it%04d.jpg',s,iter)),'jpg','Quality',90);
		if non_uniform
			plot_nonuni_kernel_write(train_ensemble_get(2,dimensions,ensemble.mx), theta_grid, fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,sprintf('s%02d_it%04d_kernel.png',s,iter)));
		else
			plot_uni_kernel_write(reshape(train_ensemble_get(2,dimensions,ensemble.mx),K,L),fullfile(OUTPUT_DIRECTORY,OUTPUT_FILENAME,sprintf('s%02d_it%04d_kernel.png',s,iter)))
		end
	end
	
    % Store current KL divergence as 'old D_val' for next iteration
    oD_val=ensemble.D_val;
    fprintf(1, [text '  Iteration %4d  Noise=%11.6e, State=%d, Repeats=%d\n'],iter,gamma_log(1,iter)^-0.5,state,tmpcount);

	% ================================================================================
    % Check if algorithm has converged
	% ================================================================================
    converged=0;
	% Algorithm has converged if change in KL divergence over last several iterations is small
    if (iter>3)
        last_dD_log=D_log(1,iter-2:iter)-D_log(1,iter-3:iter-1);
        if (all(last_dD_log<converge_criteria/1e4 & last_dD_log>-converge_criteria))
            converged=1;
        end
    end
    if (converged)
		% ================================================================================
		% If it has converged, what we do next depends on the what we were doing when it converged
		% ================================================================================
        if (state==3)
            % Have finished so exit
            break;
        elseif (state==2 & ensemble.opt_ba_sigma_2<1.1*ensemble.ba_sigma_2)
			% If we were in a state=2 iteration, and the noise parameter (inverse variance) has also stopped growing, then we are done optimising, and should just do one last iteration with state=3
			% We can't finish after a state=1 iteration, because some of the data may be randomised (see below)
            % Have converged so update noise and exit after next iteration
            state=3;
        else
			% If we weren't in a state=2 iteration, or the noise precision is still increasing, then go back for more optimisation
            % Swap between 1 and 2
            state=3-state;
        end
        last_change_iter=iter;

		% ================================================================================
        % In state=1, run through each data type and switch back on any channels of each data type of the latent variables which have converged to their prior (when q(x), which approximates the posterior, converges to the prior, it indicates the observations are not being used at all, so kick-starting q again may give it another shot). Miskin calls these "switched off components". 
		% ================================================================================
        if (state==1)
            % Set the noise directly to the optimum
            ensemble.ba_sigma_2=ensemble.opt_ba_sigma_2;
            % Randomise values for any channels that have been switched off
            % (they can always switch back off again later)
			% For each data type...
            for c=1:size(dimensions,1)
				% Get current expectations
                cmx=train_ensemble_get(c,dimensions,ensemble.mx);
                cmx2=train_ensemble_get(c,dimensions,ensemble.mx2);
				% Only randomise the channel if the options specify to restart_switched_off, if the prior is locked for this data type, and if the prior is not discrete.
				% Randomising a channel whose prior was also being optimised would probably not be a good idea, hence the 2nd condition.
                if (restart_switched_off & priorlock(c) & priortypes(c)<3)
                    %Find variance scale for this class, this is the ratio of <x>^2/<x^2>
                    %For non-rectified distributions, this ratio tends to zero as the
                    %component is switched off and the posterior tends to the prior.
                    %For rectified distributions, this ratio tends to 0.6366 as the
                    %distributions tend to rectified gaussians.
                    %Therefore randomise the components if the ratio goes below 0.7
                    cx1=train_ensemble_get(c,dimensions,ensemble.x1);
                    cx2=train_ensemble_get(c,dimensions,ensemble.x2);
                    scales=mean(cmx.^2,2)./mean(cmx2,2);
					% For each channel in this data type
                    for k=1:numchannels(c)
                        if (scales(k)<0.7)
                            %Scale of this channel is really low or C_KL is low, so randomise
                            disp(['  Reinitialising class=' int2str(c) ' k=' int2str(k)])
                            if (priortypes(c)==0)
                                %Not rectified
                                cx1(k,:)=1e4*randn(1,numelements(c)).*ceil(rand(1,numelements(c))*2);
                            else
                                %Rectified
                                cx1(k,:)=1e4*abs(randn(1,numelements(c)));
                            end
                            cx2(k,:)=1e4;
                        end
                    end
					% Put the randomised values back into the ensemble
                    ensemble.x1=train_ensemble_put(c,dimensions,ensemble.x1,cx1);
                    ensemble.x2=train_ensemble_put(c,dimensions,ensemble.x2,cx2);
                end
            end
        end
    end
end

% ================================================================================
% Shrink the logs to use the smallest possible size
% ================================================================================
D_log=D_log(:,1:iter);
gamma_log=gamma_log(1,1:iter);











