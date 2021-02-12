function [M,n] = normr(M)
    % Highly vectorized normalization of rows
    n = sqrt(sum(M.^2,2));
    M = bsxfun(@rdivide,M,n);
end