function [X_hat,P,Q,sig3] = ukf_augment(X_hat,P,Q,sig3, ii, X_hat_new,p_lmk_unc,q_lmk)
    N_prev = size(X_hat,1)-6;
    N = size(X_hat_new,1); % Number of new states
    L = size(X_hat,2);
    
    % Calcualte the segmented indices for the original matrices:
    start1 = 7;
    end1   = start1-1 + N_prev/3;
    start2 = end1+1;
    end2   = start2-1 + N_prev/3;
    start3 = end2+1;
    end3   = start3-1 + N_prev/3;
    inds1  = start1:end1;
    inds2  = start2:end2;
    inds3  = start3:end3;
    
    % Calculate the segmented indices for the new states:
    start1 = 1;
    end1   = start1-1+N/3;
    start2 = end1+1;
    end2   = start2-1+N/3;
    start3 = end2+1;
    end3   = start3-1+N/3;
    inds1_new = start1:end1;
    inds2_new = start2:end2;
    inds3_new = start3:end3;
    
    % Calculate the segmented indices for the augmented matrices:
    start1 = 7;
    end1   = start1-1 + N_prev/3;
    start2 = end1 + 1 + N/3;
    end2   = start2-1 + N_prev/3;
    start3 = end2 + 1 + N/3;
    end3   = start3-1 + N_prev/3;
    inds1_aug = start1:end1;
    inds3_aug = start2:end2;
    inds5_aug = start3:end3;
    inds2_aug = end1+1:start2-1;
    inds4_aug = end2+1:start3-1;
    inds6_aug = end3+1:7+N_prev+N-1;
    
    % Augment the state vector:
    X_hat_aug = zeros(6+N_prev+N,L);
    X_hat_aug(1:6,:)        = X_hat(1:6,:);
    X_hat_aug(inds1_aug,:)  = X_hat(inds1,:);
    X_hat_aug(inds2_aug,ii) = X_hat_new(inds1_new);
    X_hat_aug(inds3_aug,:)  = X_hat(inds2,:);
    X_hat_aug(inds4_aug,ii) = X_hat_new(inds2_new);
    X_hat_aug(inds5_aug,:)  = X_hat(inds3,:);
    X_hat_aug(inds6_aug,ii) = X_hat_new(inds3_new);
    
    % Augment the Estimation Covariance Matrix:
    SCALE = 1;
    ADD   = 5;
    P_diag = diag(P);
    P_diag_new = repelem(p_lmk_unc,N/3);
    P_diag_aug = zeros(6+N_prev+N,1);
    P_diag_aug(1:6)       = SCALE*P_diag(1:6,:)+ADD;
    P_diag_aug(inds1_aug) = P_diag(inds1,:);
    P_diag_aug(inds2_aug) = SCALE*P_diag_new(inds1_new)+ADD;
    P_diag_aug(inds3_aug) = P_diag(inds2,:);
    P_diag_aug(inds4_aug) = SCALE*P_diag_new(inds2_new)+ADD;
    P_diag_aug(inds5_aug) = P_diag(inds3,:);
    P_diag_aug(inds6_aug) = SCALE*P_diag_new(inds3_new)+ADD;
    P_aug = diag(P_diag_aug);

    % Augment the 3-sigma bound tracker:
    sig3_aug = zeros(6+N_prev+N,L);
    sig3_aug(1:6,:)        = sig3(1:6,:);
    sig3_aug(inds1_aug,:)  = sig3(inds1,:);
    sig3_aug(inds2_aug,ii) = 3*sqrt(P_diag_aug(inds1_new));
    sig3_aug(inds3_aug,:)  = sig3(inds2,:);
    sig3_aug(inds4_aug,ii) = 3*sqrt(P_diag_aug(inds2_new));
    sig3_aug(inds5_aug,:)  = sig3(inds3,:);
    sig3_aug(inds6_aug,ii) = 3*sqrt(P_diag_aug(inds3_new));
    
    % Augment the Process Noise Covariance Matrix:
    Q_diag = diag(P);
    Q_diag_new = repelem(q_lmk,N/3);
    Q_diag_aug = zeros(6+N_prev+N,1);
    Q_diag_aug(1:6)       = Q_diag(1:6,:);
    Q_diag_aug(inds1_aug) = Q_diag(inds1,:);
    Q_diag_aug(inds2_aug) = Q_diag_new(inds1_new);
    Q_diag_aug(inds3_aug) = Q_diag(inds2,:);
    Q_diag_aug(inds4_aug) = Q_diag_new(inds2_new);
    Q_diag_aug(inds5_aug) = Q_diag(inds3,:);
    Q_diag_aug(inds6_aug) = Q_diag_new(inds3_new);
    Q_aug = diag(Q_diag_aug);
        
    % Add the new initial estimate:
    X_hat = X_hat_aug;
    sig3  = sig3_aug;
    P     = P_aug;
    Q     = Q_aug;
end