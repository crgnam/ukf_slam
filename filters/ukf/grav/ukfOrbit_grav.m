function [X] = ukfOrbit_grav(dt, X, body)
    % Propagate orbit:
    X(1:6) = rk4(@gravFieldDynamics, dt, X(1:6), body);
    
    % Rotate lmks by asteroid spin rate:
    new_lmks = body.rotate(dt)'*reshape(X(7:end),[],3)';
    X(7:end) = [new_lmks(1,:)'; new_lmks(2,:)'; new_lmks(3,:)'];
end