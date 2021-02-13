classdef Camera < handle
    properties
        K
        fov
        sensor
        resolution
    end
    
    %% Constructor:
    methods (Access = public)
        function [self] = Camera(K,fov,sensor,resolution)
            self.K = K;
            self.fov = fov;
            self.sensor = sensor;
            self.resolution = resolution;
        end
    end
    
    %% Public Methods:
    methods (Access = public)
        % Function to project world points into the image space:
        function [imagePoints,inFOV] = worldToImage(self,worldPoints,rotMat,position,K2,fov2)
            % worldPoints = 3xN world points
            % rotMat   = 3x3 rotation matrix
            % position = 3x1 position vector
            % imagePoints = 3xN image points (3rd value is focal length)

            % Transform world points into camera frame:
            cameraPoints = rotMat*(worldPoints - position);

            % Calculate homogeneous coordinates:
            if nargin > 4
                homogeneous = K2*[cameraPoints; ones(1,size(cameraPoints,2))];
            else
                homogeneous = self.K*[cameraPoints; ones(1,size(cameraPoints,2))];
            end

            % Normalize the points into focal length coordinates
            imagePoints(1,:) = homogeneous(1,:)./homogeneous(3,:);
            imagePoints(2,:) = homogeneous(2,:)./homogeneous(3,:);

            % Reject points which do not fall on the sensor, or which are behind
            % the camera
            if nargin > 4
                remove1 = abs(imagePoints(1,:))>fov2(1);
                remove2 = abs(imagePoints(2,:))>fov2(2);
            else
                remove1 = abs(imagePoints(1,:))>self.fov(1);
                remove2 = abs(imagePoints(2,:))>self.fov(2);
            end
            remove3 = homogeneous(3,:)>0;
            remove = (remove1 | remove2 | remove3);

            % Pass out indices of the visible features:
            inFOV = ~remove;
        end
        
        % Function to generate rays which can be traced out from the camera
        function [rays] = generateRays(imagePoints,rotMat,K)
            % Function generates rays (as unit vectors of shape 3xN) given a
            % set of image points (2xN) and a camera projection matrix K
            raysCamFrame = [imagePoints; -K(1,1)*ones(1,size(imagePoints,2))];
            rays = normc(rotMat'*raysCamFrame);
        end
    end
end