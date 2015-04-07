function K = estimate_kernel_multi_scale(B,ks,params)
ks = ks+(1-mod(ks,2));
B = im2double(B);
%% prepare_multi_scale
scale_factor = params.scale_factor;
ks_min = params.ks_min;
[nLevels,mbsz,mksz,mlsz] = prepare_multi_scale(size(B),ks,ks_min,scale_factor);
%% initialization
Li = padarray(B,floor(ks/2),'replicate');
Ki = kernel_initialize(mksz(1,:),params.initializeMethod);
%% Begin AM
interpMethod = params.interpMethod;
nloops = params.nloops;%AM each level
eta = 1;
% begin each level
for iLevel = 1:nLevels    
    fprintf('\nLevel %d/%d ...\n',iLevel,nLevels);    
    Bi = imresize(B,mbsz(iLevel,:),interpMethod);
    Ki = imresize(Ki,mksz(iLevel,:),interpMethod);
    Ki = kernel_normalize(Ki,params.kernel_threshold); 
    Li = imresize(Li,mlsz(iLevel,:),interpMethod);         
%     figure(1),subplot(2,2,1);
%     imshow(mat2gray(imresize(K_gt,size(Ki),interpMethod)));title('Ground Truth');
%     figure(1),subplot(2,2,3);
%     imshow(mat2gray(imresize(L_gt,size(Li),interpMethod)));title('Sharp Image');
    lambda = min(params.lambda_max,params.lambda_step^(nLevels-iLevel)*params.lambda);
    Li_old = Li;    
    for jloop = 1:nloops 
        eta_L = 1/2;
        eta_K = 1/2;
        % update images
        [Li,eta_L] = update_image_armijo(Bi,Ki,Li,eta,lambda,eta_L);
        % update blur kernel
        [Ki,eta_K] = update_kernel_armijo(Bi,Ki,Li,params.gamma,eta_K);
        %normalize kernel
        Ki = kernel_normalize(Ki,params.kernel_threshold);
        if mod(jloop,100) == 0
            fprintf('%d,',jloop); 
            figure(1),subplot(2,2,2),imshow(mat2gray(Ki));
            title(sprintf('%d:%d',iLevel,jloop)); 
            figure(1),subplot(2,2,4);
            imshow(mat2gray(Li));title('Current Image');
            pause(0.05); 
        end
    end% end loop   
    Li = non_blind_decov(Bi,Ki,lambda,nloops,Li_old);
end%iLevel
K = Ki;

end