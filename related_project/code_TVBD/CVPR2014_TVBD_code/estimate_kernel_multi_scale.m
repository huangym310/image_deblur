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
        % update images
        for iter = 1:params.iters
            dLi = conv2(conv2(Li,Ki,'valid')-Bi,rot90(Ki,2),'full')+lambda*gradTVcc(Li);
            eta_Li = 5e-3*max(Li(:))/max(1e-31,max(abs(dLi(:))));
            Li = Li - eta_Li*dLi;
        end%iter
        % update blur kernel
        for iter = 1:params.iters
            dKi = conv2(rot90(Li,2),conv2(Li,Ki,'valid')-Bi,'valid');
            eta_K = 1e-3*max(Ki(:))/max(1e-31,max(abs(dKi(:))));
            Ki = Ki - eta_K*dKi;
        end%
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