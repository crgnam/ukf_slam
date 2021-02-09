function [imagePoints] = camera(rotMat,position,K,fov, worldPoints)
    % rotMat   = 3x3 rotation matrix
    % position = 3x1 position vector
    % K        = 3x3 camera projection matrix
    % fov      = 2x1 field of view limits
    % worldPoints = 3xN world points
    % imagePoints = 3xN image points (3rd value is focal length)
    
    imagePoints = K*[rotMat, position]*[worldPoints; ones(1,size(worldPoints,2))];
    imagePoints(:,abs(imagePoints(1,:))>fov(1)) = [];
    imagePoints(:,abs(imagePoints(2,:))>fov(2)) = [];
end