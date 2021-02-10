function [rays] = generateRays(imagePoints,rotMat,K)
    % Function generates rays (as unit vectors of shape 3xN) given a
    % set of image points (2xN) and a camera projection matrix K
    raysCamFrame = [imagePoints; -K(1,1)*ones(1,size(imagePoints,2))];
    rays = normc(rotMat'*raysCamFrame);
end