function [imagePoints,varargout] = camera(rotMat,position,K,fov, worldPoints)
    % rotMat   = 3x3 rotation matrix
    % position = 3x1 position vector
    % K        = 3x3 camera projection matrix
    % fov      = 2x1 field of view limits
    % worldPoints = 3xN world points
    % imagePoints = 3xN image points (3rd value is focal length)
    
    % Calculate the homogeneous coordinates:
    cameraPoints = rotMat*(worldPoints - position);
    homogeneous = K*[cameraPoints; ones(1,size(cameraPoints,2))];
    
    % Normalize the points into focal length coordinates
    imagePoints(1,:) = homogeneous(1,:)./homogeneous(3,:);
    imagePoints(2,:) = homogeneous(2,:)./homogeneous(3,:);
    
    % Reject points which do not fall on the sensor, or which are behind
    % the camera
    inds1 = abs(imagePoints(1,:))>fov(1);
    inds2 = abs(imagePoints(2,:))>fov(2);
    inds3 = homogeneous(3,:)>0;
    inds = (inds1 | inds2 | inds3);
    imagePoints(:,inds) = [];

    % Pass out keep indices for color rendering (if requested)
    if nargout == 2
        varargout{1} = ~inds;
    end
end