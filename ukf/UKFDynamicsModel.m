function [X] = UKFDynamicsModel(dt, X, varargin)
    X(1:6) = rk4(@orbitalDynamics, dt, X(1:6), varargin{:});
end