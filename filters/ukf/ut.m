%% Unscented Transform:
% This is a simple implementation of the Unscented Transform (UT).  It is
% primarily meant for use in propagating states/covariances in the UKF.
%
% Inputs:
%   systemModel - function handle to the system model used for propagation
%   dt     - time step size
%   sigmas - sigma points to be propagated
%   Wm     - mean weights
%   Wc     - covariance weights
%   R      - additional covariance
%   n      - number of sigma points
%   varargin - arguments needed for the system model
%
% Outputs:
%   mu - updated mean
%   P  - a posteriori distribution
%   deviations - (mu - x_hat)
%   sigmas     - propagated sigma points
%
% Requirements:
%     aposteriori_distribution() - constructs the aposteriori distribution
%
% For reference, see:
%    Optimal Estimation of Dynamic Systems - Crassidis
%
% Chris Gnam - 2019

function [mu, P, deviations, sigmas_out] = ut(systemModel, dt, sigmas, Wm, Wc, R, n_out, varargin)
    
    num_sigmas = size(sigmas,2);
    sigmas_out = zeros(n_out, num_sigmas);
    
    % Propagate sigma points through dynamics:
    for ii = 1:num_sigmas
        sigmas_out(:,ii) = systemModel(dt, sigmas(:,ii), varargin{:});
    end
    
    % Calculate new mean:
    mu = sigmas_out*Wm'; 
    
    % Recover a posteriori distribution:
    [P, deviations] = aposteriori_distribution(sigmas_out, mu, num_sigmas, Wc, R);
end