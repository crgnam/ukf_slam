function [X] = ukfOrbit_bennu_rotation(dt, X, varargin)
    X(1:6) = rk4(@orbitalDynamics, dt, X(1:6), varargin{:});
    X(9) = X(9) + X(10)*dt;
end