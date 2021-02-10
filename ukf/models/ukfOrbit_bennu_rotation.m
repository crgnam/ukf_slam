function [X] = ukfOrbit_bennu_rotation(dt, X, varargin)
    % Propagate the orbit forward:
    X(1:6) = rk4(@orbitalDynamics, dt, X(1:6), varargin{:});
    
    % Rotate bennu:    
    X(9) = X(9) + X(10)*dt;
    bennu_rotMat = ea2rotmat123(X(7),X(8),X(10)*dt);
    new_lmks = bennu_rotMat*reshape(X(11:end),[],3)';
    X(11:end) = [new_lmks(1,:)'; new_lmks(2,:)'; new_lmks(3,:)'];
end