function [predicted_measurement] = ukfCamera(~,X,rotMat,K, avails)
    r = X(1:3);
%     v = X(4:6);
%     N = (length(X)-6)/3;
    imagePoints = camera(rotMat,r,K,[nan nan], reshape(X(7:end),[],3)');
    predicted_measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    predicted_measurement = predicted_measurement(avails);
end