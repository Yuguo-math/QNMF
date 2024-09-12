function [AA_hat E_hat iter] = inexact_alm_NMF_MC(D,Support, C, myeps, tol, maxIter)


UpdatingSupport = ~Support;

CUpdatingSupport = UpdatingSupport(:,:,1);
CUpdatingSupport = [CUpdatingSupport;CUpdatingSupport;CUpdatingSupport;CUpdatingSupport];


if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
CD = q2cplx( double2q(D) );
[m n] = size(CD);

Y = CD;
norm_two = max(svd(Y, 'econ'));
norm_inf = norm( Y(:), inf) ;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);


mu = 1/norm_two;  
mu_bar = mu * 1e7;
 
rho = 1.05;         
d_norm = norm(CD, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;

while ~converged       
    iter = iter + 1;
    
    
    temp_T = CD - A_hat + (1/mu)*Y;    
    E_hat = CUpdatingSupport.*temp_T;  
    
    [U S V] = svd(CD - E_hat + (1/mu)*Y, 'econ');

    
    alpha = 4;


    [tempDiagS,svp] = ClosedQNMF(S, alpha, C/mu);

    A_hat = U(:,1:svp)*diag(tempDiagS)*V(:,1:svp)';  % X

    %%%%%%%%%%%%
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end

    %%%%%% 
    total_svd = total_svd + 1;
    
    
    Z = CD - A_hat - E_hat;
    
    
    Y = Y + mu*Z;

    mu = min(mu*rho, mu_bar);


    AA_hat = double2q(q2cplx(A_hat, 'inverse'), 'inverse');
    AA_hat = real(AA_hat);
    imshow(AA_hat/255)
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    
    if ~converged && iter >= maxIter
        converged = 1 ;       
    end
end
