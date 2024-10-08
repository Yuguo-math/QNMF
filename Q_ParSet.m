function  [par]=Q_ParSet(nSig)

par.nSig      =   nSig;                                   % Variance of the noise image
par.SearchWin =   30;                                     % Non-local patch searching window
par.delta     =   0.1;                                    % Parameter between each iter
par.c         =   2.5*sqrt(2); %2.8 1.3*sqrt(2);%%sqrt(2) % Constant num for the weight vector
par.Innerloop =   2;                                    % InnerLoop Num of between re-blockmatching
par.ReWeiIter =   1;%%3

if nSig<=5
    par.patsize       =   4;                            % Patch size
    par.patnum        =   80;                           % Initial Non-local Patch number
    par.Iter          =   4;                            % total iter numbers
    par.lamada        =   0.54;                         % Noise estimete parameter
    
elseif nSig<=20
    par.patsize       =   4;                            % Patch size
    par.patnum        =   80;                           % Initial Non-local Patch number
    par.Iter          =   6;                            % total iter numbers
    par.lamada        =   0.54;                         % Noise estimete parameter

elseif nSig <= 40
    par.patsize       =   5;  
    par.patnum        =   90;
    par.Iter          =   8;
    par.lamada        =   0.56; 

elseif nSig<=60
    par.patsize       =   5; 
    par.patnum        =   120;
    par.Iter          =   9;  
    par.lamada        =   0.58; 
    
else
    par.patsize       =   5;  
    par.patnum        =   140;
    par.Iter          =   10;
    par.lamada        =   0.58; 
end

par.step      =   floor((par.patsize)/2-1);         
