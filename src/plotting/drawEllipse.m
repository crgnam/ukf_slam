function [varargout] = drawEllipse(x0,y0,a,b,R,varargin)
    % Calculate ellipse points:
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = x0 + a*cos(theta_r);
    ellipse_y_r     = y0 + b*sin(theta_r);
    rotated_ellipse = R*[ellipse_x_r;ellipse_y_r];
    
    % Plot the ellipse:
    h = plot(rotated_ellipse(1,:),rotated_ellipse(2,:),'color',varargin{:});
    if nargout == 1
        varargout = {h};
    end
end