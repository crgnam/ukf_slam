function [predicted_measurement] = ukfCamera_bennu_rotation(~,X,rotMat,K, avails)
    r = X(1:3);
%     v = X(4:6);
    theta = X(7);
    phi = X(8);
    psi = X(9);
%     w = X(10);
    
    bennu_rotMat = ea2rotmat123(theta,phi,psi);
    
    imagePoints = camera(rotMat,r,K,[nan nan], bennu_rotMat*reshape(X(11:end),[],3)');
    predicted_measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    predicted_measurement = predicted_measurement(avails);
end