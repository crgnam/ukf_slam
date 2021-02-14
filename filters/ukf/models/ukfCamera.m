function [predicted_measurement] = ukfCamera(~,X,satellite,avails)
    r = X(1:3);
%     v = X(4:6);
    [imagePoints,~] = satellite.camera.worldToImage(reshape(X(7:end),[],3)',...
                                                    satellite.rotmat,r);
    predicted_measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    predicted_measurement = predicted_measurement(avails);
end