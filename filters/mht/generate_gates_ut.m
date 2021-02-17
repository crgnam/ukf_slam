function [gates] = generate_gates_ut(X_hat,P,Q,R, dt, ukf_args, measAvails)   
    % Extract UKF Information:
    dynamics_args = ukf_args{1};
    measurement_args = ukf_args{2};
    
    % Generate sigma points:
    [SIGMAS,Wm,Wc,L] = u_sigmas(X_hat,P,1,1,1);
    [~,~,~,SIGMAS] = ut(@ukfOrbit, dt, SIGMAS, Wm, Wc, Q, L, dynamics_args{:});
    num_sigmas = size(SIGMAS,2);
    nn = sum(measAvails,1);
    sigmas_out = zeros(nn, num_sigmas);
    
    % Formulate measurement covariance matrix for available measurements:
    Rdiag = diag(R);
    Rdiag = Rdiag(measAvails);
    R = diag(Rdiag);   
    
    % Projection sigma points into the image space:
    for ii = 1:num_sigmas
        sigmas_out(:,ii) = ukfCamera(0, SIGMAS(:,ii), measurement_args{:}, measAvails);
    end
    
    % Calculate new mean:
    mu = sigmas_out*Wm'; 
    
    % Recover a posteriori distribution:
    [P,~] = aposteriori_distribution(sigmas_out, mu, num_sigmas, Wc, R);
    
    % Format gates:
    P_diag = diag(P);
    mux = mu(1:nn/2);
    muy = mu(1+nn/2:nn);
    sig3x = 3*sqrt(P_diag(1:nn/2));
    sig3y = 3*sqrt(P_diag(1+nn/2:nn));
    gates = [mux,muy,sig3x,sig3y];
end