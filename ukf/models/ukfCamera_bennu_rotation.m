function [predicted_measurement] = ukfCamera_bennu_rotation(~,X,rotMat,K, avails)
    r = X(1:3);
%     v = X(4:6);
%     N = (length(X)-6)/3;
    imagePoints = camera(rotMat,r,K,[nan nan], reshape(X(11:end),[],3)',nan,0);
    predicted_measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    predicted_measurement = predicted_measurement(avails);
end