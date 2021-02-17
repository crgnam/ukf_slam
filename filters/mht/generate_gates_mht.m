function [gates] = generate_gates_mht(X_hat,P)   
    % Extract predicted measurement data:
    meas = X_hat(7:end);
    P_meas = P(7:end,7:end);
    nn = size(meas,1)/2;

    % This code came directly from Scott's MHT code:
    p_D = 1/nn; % Probability that a detection is made if a true target is present
    B = 1; %?????
    M = 2; % Dimension of measurements
        
    % Calculate the gate radius for each measurement:
    G = zeros(nn,1);
    for ii = 1:nn
        G(ii) = 2*log(p_D/((1-p_D)*(2*pi)^(M/2)*B*sqrt(det( P_meas((2*ii:2*ii+1)-1,(2*ii:2*ii+1)-1)) )));
    end
    
    % Format gates:
    mux = meas(1:nn);
    muy = meas(1+nn:end);
    gates = [mux,muy,G,G];
end