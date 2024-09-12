
function  [X] =  QNMF( Y, C, NSig, m, Iter)

    % 转四元数后cplx ！
    cY = q2cplx( double2q(Y) );

    [U,SigmaY,V] =   svd(full(cY),'econ');    

    PatNum  = size(Y,2);
    mNSig   = sqrt(sum(NSig(:).^ 2)/length(NSig(:)));
    TempC   = 2*C*sqrt(PatNum)*mNSig^2;
    
    alpha = 4;
    rho = 1;

    [SigmaX,svp] = ClosedQNMF(SigmaY, alpha, TempC/rho);

    
    [~, h] = size(U); 
    if svp>h        
        svp = h;
        SigmaX = SigmaX(1:h,:);
    end
    
    cX =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     

    X= double2q(q2cplx(cX, 'inverse'), 'inverse');
    X= real(X) + m;
return;
