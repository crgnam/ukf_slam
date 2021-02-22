function [dX] = gravFieldDynamics(~,X,body,varargin)
    a = body.gravityField.acceleration(X(1:3), body.inert2body, varargin{:});
    dX = [X(4:6); a];
end