function [imagePoints,varargout] = camera(rotMat,position,K,fov, worldPoints,normals, sig_meas)
    % rotMat   = 3x3 rotation matrix
    % position = 3x1 position vector
    % K        = 3x3 camera projection matrix
    % fov      = 2x1 field of view limits
    % worldPoints = 3xN world points
    % sig_pix  = Standard deviation of pixel noise
    % imagePoints = 3xN image points (3rd value is focal length)
    
    % Transform world points into camera frame:
    cameraPoints = rotMat*(worldPoints - position);
    
    % Remove features pointed away from the camera:
    if ~any(isnan(normals)) % normals == nan allows us to bypass this step
        cameraNorms  = rotMat*normals;
        inds0 = cameraNorms(3,:)<0;
    else
        inds0 = false(1,size(cameraPoints,2));
    end

    % Calculate homogeneous coordinates:
    homogeneous = K*[cameraPoints; ones(1,size(cameraPoints,2))];
    
    % Normalize the points into focal length coordinates
    imagePoints(1,:) = homogeneous(1,:)./homogeneous(3,:);
    imagePoints(2,:) = homogeneous(2,:)./homogeneous(3,:);
    
    % Reject points which do not fall on the sensor, or which are behind
    % the camera
    if ~any(isnan(fov))
        inds1 = abs(imagePoints(1,:))>fov(1);
        inds2 = abs(imagePoints(2,:))>fov(2);
        inds3 = homogeneous(3,:)>0;
        inds = (inds0 | inds1 | inds2 | inds3);
        imagePoints(:,inds) = [];
    else
        inds = nan;
    end
    
    % Flip the Image along both axes:
    imagePoints = -imagePoints;
    imagePoints = imagePoints + sig_meas*randn(size(imagePoints));

    % Pass out keep indices for color rendering (if requested)
    if nargout == 2
        varargout{1} = ~inds;
    end
end