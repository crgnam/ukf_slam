%% Aposteriori Distribution:
% Inputs:
%   sigmas - sigma points
%   mu     - a posteriori updated mean
%   n      - number of sigma points
%   Wc     - covariance weights
%   R      - additional covariance
%
% Outputs:
%   P  - a posteriori distribution
%   deviations - (mu - x_hat)
%
% For reference, see:
%    Optimal Estimation of Dynamic Systems - Crassidis
%
% Chris Gnam - 2019

function [P, deviations] = aposteriori_distribution(sigmas, mu, n, Wc, R)
    % Calculate the deviations:
    deviations = sigmas - mu(:,ones(1,n));
    
    % Recover distribution:
    P  = deviations*diag(Wc)*deviations' + R; 
end