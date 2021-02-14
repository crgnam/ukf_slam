function [h] = drawSphere(xc,yc,zc,R,varargin)
    [X,Y,Z] = sphere;
    
    % Scale to desire radius.
    X = X*R;
    Y = Y*R;
    Z = Z*R;
    
    % Plot as surface.
    h = surf(X+xc, Y+yc, Z+zc,varargin{:});
end