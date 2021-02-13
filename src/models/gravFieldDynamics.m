function [dX] = gravFieldDynamics(~,X,body)
    a = body.gravityField.acceleration(X(1:3), body.inert2body);
    dX = [X(4:6); a];
end