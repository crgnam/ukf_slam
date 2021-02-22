function [X] = ukfOrbit_grav(dt, X, body, num_uwb)
    % Extract vehicle states:
    rv = X(1:6);
    
    % Extract the states of the uwb transceivers:
    r_uwb = X(7:6+3*num_uwb);
    v_uwb = X(7+3*num_uwb:6+6*num_uwb);
    
    % Extract the gravity field information:
    mu = X(7+6*num_uwb);
    coeffs_vec = X(8+6*num_uwb:end);
    [Cnm,Snm] = vec2coeffs(coeffs_vec);
    
    % Propagate all orbits:
    rv = rk4(@gravFieldDynamics, dt, rv, body, Cnm,Snm,mu);
    r_uwb = reshape(r_uwb,3,[]);
    v_uwb = reshape(v_uwb,3,[]);
    for ii = 1:num_uwb
        rv_uwb = rk4(@gravFieldDynamics,dt,[r_uwb(:,ii); v_uwb(:,ii)],body, Cnm,Snm,mu);
        r_uwb(:,ii) = rv_uwb(1:3);
        v_uwb(:,ii) = rv_uwb(4:6);
    end
    
    % Store outputs:
    X = [rv;r_uwb(:);v_uwb(:); mu; coeffs_vec];
end