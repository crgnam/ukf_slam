function [dX] = orbitalDynamics(~, X, mu)
    % State space representation of two-body orbital dynamics:
    r = X(1:3);
    v = X(4:6);
    dX = [v; -(mu*r)/norm(r)^3];
end