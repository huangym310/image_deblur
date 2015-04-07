function [fp, Mp, Np, MKp, NKp, lambdas, scales] = buildPyramid(f, ...
                                                      MK, NK, lambda, ...
                                                      lambdaMultiplier, ...
                                                      interpMethod, ...
                                                      scaleMultiplier, ...
                                                      largestLambda)
% Author: Daniele Perrone, perrone@iam.unibe.ch
% Copyright (C) 2014 Daniele Perrone.  All rights reserved.

    [M, N, ~] = size(f);
    
    smallestScale = 3;
    
    if (~exist('largestLambda','var'))
        largestLambda = 1.1e-1;
    end
    
    if (~exist('scaleMultiplier','var'))
        scaleMultiplier = 1.1;
    end
    
    scales = 1;
  
    
    fp{1,1}  = f;
    Mp{1,1} = M;
    Np{1,1} = N;
    MKp{1,1} = MK;
    NKp{1,1} = NK;
    
    lambdas(1) = lambda;
    
    while (MKp{scales,1} > smallestScale && NKp{scales,1} > smallestScale...
            && lambdas(scales) * lambdaMultiplier < largestLambda)
        scales = scales + 1;
        
        % Compute lambda value for the current scale
        lambdas(scales) = lambdas(scales - 1) * lambdaMultiplier;
       
        MKp{scales,1} = round(MKp{scales - 1,1} / scaleMultiplier);
        NKp{scales,1} = round(NKp{scales - 1,1} / scaleMultiplier);
        
        
        % Makes kernel dimension odd
        if (mod(MKp{scales,1},2) == 0)
            MKp{scales,1} = MKp{scales,1} - 1;
        end
        if (mod(NKp{scales,1},2) == 0)
            NKp{scales,1} = NKp{scales,1} - 1;
        end
        
        if (NKp{scales,1} == NKp{scales-1,1})
            NKp{scales,1} = NKp{scales,1} - 2;
        end
        if (MKp{scales,1} == MKp{scales-1,1})
            MKp{scales,1} = MKp{scales,1} - 2;
        end
        
        if (NKp{scales,1} < smallestScale)
            NKp{scales,1} = smallestScale;
        end
        
        if (MKp{scales,1} < smallestScale)
            MKp{scales,1} = smallestScale;
        end
        
        % Correct scaleFactor for kernel dimension correction
        factorM = MKp{scales - 1,1}/MKp{scales,1};
        factorN = NKp{scales - 1,1}/NKp{scales,1};
        
        Mp{scales,1} = round(Mp{scales - 1,1} / factorM);
        Np{scales,1} = round(Np{scales - 1,1} / factorN);
        
        
        % Makes image dimension odd
        if (mod(Mp{scales,1},2) == 0)
            Mp{scales,1} = Mp{scales,1} - 1;
        end
        if (mod(Np{scales,1},2) == 0)
            Np{scales,1} = Np{scales,1} - 1;
        end
        

        fp{scales,1} = imresize(f,[Mp{scales,1} Np{scales,1}], 'Method', interpMethod);       
    end
end