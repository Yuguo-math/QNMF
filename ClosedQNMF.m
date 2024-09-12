function [SigmaZ, svp] = ClosedQNMF(SigmaY, alpha, C)
    
    SigmaY = diag(SigmaY);

    norm_inf = SigmaY(1);
    Lambda = 2*sqrt(C);
    
    if (norm_inf>Lambda)

        z = max(abs(SigmaY)-Lambda, 0);
        svp = length(find(z > 0));
        z = z(1:svp) + Lambda/2;
    
        norm_2 = norm(z, 2);
        K = 1 + alpha*( Lambda/2 )/norm_2;
        SigmaZ = K*z ;


        over = find(SigmaZ > SigmaY(1:svp)); 
        SigmaZ(over) = SigmaY(over);

        ind = find(z > 0); 
        svp = length(ind);
        SigmaZ = SigmaZ(ind);
        
    else 
        SigmaZ = [];    svp = 0;
    end

return;