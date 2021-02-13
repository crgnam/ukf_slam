classdef Spacecraft < handle
    properties
        % States of the spacecraft:
        r
        v
        rotmat
        
        % camera object:
        camera
        
        % visualization stuff:
        plotted_spacecraft = false;
        xyz
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = Spacecraft(r,v,rotmat,camera)
            self.r = r;
            self.v = v;
            self.rotmat = rotmat;
            self.camera = camera;
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Propagate the spacecraft forward in time:
        function [self] = propagate(self,dt,body)
            X = rk4(@gravFieldDynamics,dt,[self.r; self.v],body);
            self.r = X(1:3);
            self.v = X(4:6);
            
            % Re-point nadir (or replace with rotational
            % dynamics/controller if needed in the future):
            self.rotmat = nadir(self.r,self.v);
        end
        
        % Take a measurement given a body with a set of lmks:
        function [imagePoints,visible] = image(self,body,sig_meas)
            % Project the body landmark locations to the image space:
            [imagePoints,inFOV] = self.camera.worldToImage(body.lmks_i,self.rotmat,self.r);
            
            % Detect points that are pointed towards the camera:
            cameraNorms = self.rotmat*body.lmk_norms_i;
            visible = cameraNorms(3,:)>0;
            
            % Detect points that are illuminated:
            if norm(body.sun_vec) == 0
                illuminated = true;
            else
                dir_sgn = sum(body.sun_vec.*body.lmk_norms_i,1);
                illuminated = dir_sgn>0;
            end
            
            % Select only points that are in FOV, illuminated, and visible:
            visible = inFOV & visible & illuminated;
            imagePoints = imagePoints(:,visible);
            
            % Add noise to the measurement image points:
            imagePoints = imagePoints + sig_meas*randn(size(imagePoints));
        end
    end
    
    %% Public Methods for Visualizations
    methods (Access = public)
        % Draw the spacecraft's current attitude:
        function [] = draw(self,varargin)
            if isnumeric(varargin{1})
                scale = varargin{1};
                if nargin > 2
                    varargin = varargin(2:end);
                end
            else
                scale = 1;
            end
            
            rotMat = scale*self.rotmat;
            if ~self.plotted_spacecraft
                self.xyz = drawOrientation(self.r,rotMat,varargin{:});
                self.plotted_spacecraft = true;
            else
                self.xyz = updateOrientation(self.xyz,self.r,rotMat);
            end
        end
    end
end