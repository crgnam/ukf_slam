classdef Spacecraft < handle
    properties
        % States of the spacecraft:
        r
        v
        rotmat
        rotmat_hat
        
        % camera object:
        camera
        
        % visualization stuff:
        plotted_spacecraft = false;
        plotted_spacecraftEstimate = false;
        plotted_traj = false;
        traj_handle
        xyz
        xyz_hat
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
        function [imagePoints,visible,lmk_inds] = imageBody(self,body,sig_meas)
            % Project the body landmark locations to the image space:
            [imagePoints,inFOV] = self.camera.worldToImage(body.lmks_i,self.rotmat,self.r);
            
            % Detect points that are pointed towards the camera:
            cameraNorms = self.rotmat*body.lmk_norms_i;
            inview = cameraNorms(3,:)>0;
            
            % Detect points which meet with viewing angle requirement:
            rays   = normc(body.lmks_i - self.r);
            angled = acos(sum(rays.*-body.lmk_norms_i,1)) < self.camera.max_va;
            
            % Detect points that are illuminated:
            if norm(body.sun_vec) == 0
                illuminated = true;
            else
                dir_sgn = sum(body.sun_vec.*body.lmk_norms_i,1);
                illuminated = dir_sgn>0;
            end
            
            % Select only points that are in FOV, illuminated, and visible:
            visible = inFOV & inview & angled & illuminated;
            lmk_inds = 1:size(body.lmks_i,2);
            lmk_inds = lmk_inds(visible);
            imagePoints = imagePoints(:,visible);
            imagePoints = -imagePoints; % Import for ray tracing steps
            
            % Add noise to the measurement image points:
            imagePoints = imagePoints + sig_meas*randn(size(imagePoints));
        end
        
        % Take a measurement given a set of 3d locations:
        function [imagePoints,visible] = image(self,points,body,sig_meas)
            % Project the body landmark locations to the image space:
            [imagePoints,inFOV] = self.camera.worldToImage(points,self.rotmat,self.r);
            
            % Detect points that are not occluded by the body:
            cameraNorms = self.rotmat*body.lmk_norms_i;
            inview = cameraNorms(3,:)>0;
            
            % Detect points which meet with viewing angle requirement:
            rays   = normc(body.lmks_i - self.r);
            angled = acos(sum(rays.*-body.lmk_norms_i,1)) < self.camera.max_va;
            
            % Detect points that are illuminated:
            if norm(body.sun_vec) == 0
                illuminated = true;
            else
                dir_sgn = sum(body.sun_vec.*body.lmk_norms_i,1);
                illuminated = dir_sgn>0;
            end
            
            % Select only points that are in FOV, illuminated, and visible:
            visible = inFOV & inview & angled & illuminated;
            lmk_inds = 1:size(body.lmks_i,2);
            lmk_inds = lmk_inds(visible);
            imagePoints = imagePoints(:,visible);
            imagePoints = -imagePoints; % Import for ray tracing steps
            
            % Add noise to the measurement image points:
            imagePoints = imagePoints + sig_meas*randn(size(imagePoints));
        end
    end
    
    %% Public Methods for Visualizations
    methods (Access = public)
        % Reset all drawing data:
        function [self] = reset(self)
            self.plotted_spacecraft = false;
            self.plotted_spacecraftEstimate = false;
        end
        
        % Draw the spacecraft's current attitude:
        function [] = draw(self,varargin)
            if nargin > 1
                if isnumeric(varargin{1})
                    scale = varargin{1};
                    if nargin > 2
                        varargin = varargin(2:end);
                    else
                        varargin = {};
                    end
                else
                    scale = 1;
                end
            else
                scale = 1;
            end
            
            rotMat = scale*self.rotmat;
            if ~self.plotted_spacecraft
                self.xyz = drawOrientation(self.r,rotMat,'-',varargin{:});
                self.plotted_spacecraft = true;
            else
                self.xyz = updateOrientation(self.xyz,self.r,rotMat);
            end
        end
        
        % Draw the trajectory history of the spacecraft:
        function [self] = drawTraj(self,r_hist,varargin)
            if ~self.plotted_traj
                self.traj_handle = plot3(r_hist(1,:),...
                                         r_hist(2,:),....
                                         r_hist(3,:),...
                                         varargin{:}); hold on
                self.plotted_traj = true;
            else
                set(self.traj_handle,'XDAta',r_hist(1,:),...
                                     'YData',r_hist(2,:),...
                                     'ZData',r_hist(3,:))
            end
        end
        
        % Draw a depiction of the spacecrat's estimated position:
        function [] = drawEstimate(self,r_hat,varargin)
            if nargin > 1
                if isnumeric(varargin{1})
                    scale = varargin{1};
                    if nargin > 2
                        varargin = varargin(2:end);
                    else
                        varargin = {};
                    end
                else
                    scale = 1;
                end
            end
            
            rotMat = scale*self.rotmat_hat;
            if ~self.plotted_spacecraftEstimate
                self.xyz_hat = drawOrientation(r_hat,rotMat,'--',varargin{:});
                self.plotted_spacecraftEstimate = true;
            else
                self.xyz_hat = updateOrientation(self.xyz_hat,r_hat,rotMat);
            end
        end
        
        % Draw the spacecraft's ray projections of image points:
        function [h] = drawRays(self,imagePoints,radius_estimate,radius_max,r_hat)
            % Calculate distance:
            distance = norm(r_hat);
    
            % Generate the rays:
            rays = self.camera.generateRays(imagePoints,self.rotmat_hat);
            origins = zeros(size(rays',1),3);
            rays = rays';
            
            % Trace the rays to find ray-sphere intersections:
            lines  = [origins rays];
            sphere_pos = self.rotmat_hat'*[0;0;-distance];
            sphere = [sphere_pos' radius_estimate];
            intersects = intersectLineSphere(lines, sphere);
                      
            if any(isnan(intersects))
%                 intersects = origins + distance*rays;
                sphere = [sphere_pos' (radius_max+radius_estimate)/2];
                intersects = intersectLineSphere(lines, sphere);
            end
            intersects = getClosestIntersects(intersects,r_hat);
            intersects = intersects - sphere_pos';
            
            % Plot the data:
            quiver3(origins(:,1)+r_hat(1),origins(:,2)+r_hat(2),origins(:,3)+r_hat(3),...
            distance*rays(:,1),distance*rays(:,2),distance*rays(:,3),'r'); hold on
            drawSphere(0,0,0,sphere(4),...
                       'FaceColor',[0 0 1],'EdgeColor','None','FaceLighting','gouraud',...
                       'FaceAlpha',0.2,'AmbientStrength',1);
            plot3(intersects(:,1),intersects(:,2),intersects(:,3),'.m','MarkerSize',30)
            axis equal
            grid on
            
        end
    end
end