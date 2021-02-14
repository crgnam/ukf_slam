function [T, p] = cholcov_codegen(Sigma)

% Test for square, symmetric
[n,m] = size(Sigma);
wassparse = issparse(Sigma);
tol = 10*eps(max(abs(diag(Sigma))));
if (n == m) && all(all(abs(Sigma - Sigma') < tol))
    [T,p] = chol(Sigma);
    
    if p > 0
        % Test for positive definiteness
        % Can get factors of the form Sigma==T'*T using the eigenvalue
        % decomposition of a symmetric matrix, so long as the matrix
        % is positive semi-definite.
        [U,D] = eig(full((Sigma+Sigma')/2));
        
        % Pick eigenvector direction so max abs coordinate is positive
        [~,maxind] = max(abs(U),[],1);
        negloc = (U(maxind + (0:n:(m-1)*n)) < 0);
        U(:,negloc) = -U(:,negloc);
        
        D = diag(D);
        tol = eps(max(D)) * length(D);
        t = (abs(D) > tol);
        D = D(t);
        p = sum(D<0); % number of negative eigenvalues
        
        if (p==0)
            T = diag(sqrt(D)) * U(:,t)';
        else
            T = zeros(0,'like',Sigma);
        end
    end
    
else
    T = zeros(0,'like',Sigma);
    p = nan('like',Sigma);
end

if wassparse
    T = sparse(T);
end
