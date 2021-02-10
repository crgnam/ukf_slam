%% Unscented Kalman Filter:
% This is a simple implementation of a general purpose Unscented Kalman
% Filter (UKF).  It is meant to be a simple tool for quick and easy
% filtering of a wide range of systems.
%
% Improvement may be made, including adding the ability to handle different
% types of measurements, or varying numbers of measurements per time step
%
% Inputs:
%   dynamics  - function handle to the system dynamics
%   measModel - function handle to the measurement model
%   X_hat     - current state estimate
%   dt        - time step size
%   P, Q, R   - vovariance, Process Noise, and Measurement Covariances
%   measAvail - is a measurement available
%   meas      - measurement vector
%   alpha, beta, k - UKF tunable parameters
%   varargin  - cell array containing dynamics_args, and measurement_args,
%               both of which are also cell arrays, and correspond to the
%               necessary inputs for the respective function handles
%
% Outputs:
%   X_hat - Updated state estimate
%   P     - Updated covariance
%
% Requirements:
%     UT() - Unscented Transform function
%
% For reference, see:
%    Optimal Estimation of Dynamic Systems - Crassidis
%
% Chris Gnam - 2019

function [X_hat, P] = ukf(dynamics, measModel, X_hat, dt,...
                          P, Q, R, measAvails, meas,...
                          alpha, beta, k, varargin)
                      
    % Process the variable inputs:
    dynamics_args = varargin{1};
    measurement_args = varargin{2};
    
    % Extract important values:
    L = size(X_hat,1);
    n = sum(measAvails,1);

    % Initial calculations:
    lambda = alpha^2*(L + k) - L;
    gamma = sqrt(L + lambda);

    % Weights (3.259):
    Wm = [lambda/(L+lambda), 1/(2*(L+lambda)) + zeros(1,2*L)]; %Mean
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - (alpha^2) + beta); %Covariance
    
    % Calculate matrix square root via Choleksy Decomposition:
    sig = gamma*chol(P)';

    % Generate set of sigma points, with estimate as mean (3.253):
    SIGMAS = [X_hat, X_hat(:, ones(1,L))+sig, X_hat(:, ones(1,L))-sig];
    
    % Propagate estimate through dynamics:
    [X_hat, P, X_dev, SIGMAS] = ut(dynamics, dt, SIGMAS, Wm, Wc, Q, L, dynamics_args{:});
    
    % Establish which measurements to process:
    if any(measAvails)
        measAvail = true;
        Rdiag = diag(R);
        Rdiag = Rdiag(measAvails);
        R = diag(Rdiag);
        meas = meas(measAvails);
        
    else
        measAvail = false;
    end
    
    % If measurement is available, perform kalman update:
    measurement_args{end+1} = measAvails;
    if measAvail
        % Calculate predicted measurement:
        [y_hat, Pyy, y_dev]  = ut(measModel, dt, SIGMAS, Wm, Wc, R, n, measurement_args{:});
        
        % Calculate cross-correlation covariance (3.266):
        Pxy = X_dev*diag(Wc)*y_dev';

        % Kalman Gain (3.251):
        K   = Pxy*Pyy^-1;
        
        % State and Covariance Update (3.249):
%         disp(meas-y_hat)
        X_hat = X_hat + K*(meas - y_hat);
        P = P - K*Pxy';
    end
    
    % Make sure covariance is positive semidefinite.
%     try
        P = posSemiDefCov(P);
%     catch
%         disp('hi')
%     end
end 