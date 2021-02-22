function [predicted_measurement] = ukfCamera_grav(~,X,satellite,estimated_map,num_uwb,avails)
    % Extract vehicle states:
    r = X(1:3);
%     v = X(4:6);
    
    % Extract the states of the uwb transceivers:
    r_uwb = reshape(X(7:6+3*num_uwb), 3,[]);
%     v_uwb = X(7+3*num_uwb:6+6*num_uwb);

    % Predict LMK measurements using estimated map:
    [lmk_imagePoints,~] = satellite.camera.worldToImage(estimated_map,...
                                                        satellite.rotmat,r);
    lmk_imagePoints = -lmk_imagePoints;
    
    % Predict UWB optical measurements:
    [uwb_imagePoints,~] = satellite.camera.worldToImage(r_uwb,...
                                                        satellite.rotmat,r);
    uwb_imagePoints = -uwb_imagePoints;
    
    % Predict UWB ranging measurements:
    vec_pairs = combvec([r, r_uwb]);
    [~,uwb_ranges] = normr(vec_pairs(:,4:6) - vec_pairs(:,1:3));
    
    % Store all as predicted_measurement:
    predicted_measurement = [lmk_imagePoints(1,:)'; lmk_imagePoints(2,:)';
                             uwb_imagePoints(1,:)'; uwb_imagePoints(2,:)';
                             uwb_ranges];
    predicted_measurement = predicted_measurement(avails);
end