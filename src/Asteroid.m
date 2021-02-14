classdef Asteroid < handle
    properties 
        % Shape model:
        faces
        verts_i
        verts_b
        
        % Landmarks:
        lmks_i % landmark inertial locations
        lmks_b % landmark body locations
        lmk_norms_i % landmark normals inertial
        lmk_norms_b % landmark normals body
        lmks_obs % boolean aray for if landmark has ever been observed
        lmks_lbl % label for landmarks in order of first observation
        
        % Environment:
        inert2body
        radius
        radius_max
        rotation_axis
        rotation_rate
        sun_vec
        
        % Gravity parameters:
        gravityField
        Cnm
        Snm
        mu
        
        % Plotting information:
        ax
        ptch
        lght
        plotted_body = false;
        plotted_lmks = false;
        lmk_vis
        lmk_inv_obs
        lmk_inv_unobs
        
        % Estimated (for plotting):
        lmks_hat
        plotted_lmks_hat = false;
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = Asteroid(obj_file, grav_file, inert2body,...
                                   rotation_rate,sun_vec, num_lmks, scale)
            % Get orientation:
            self.inert2body = inert2body;
            
            % Read in shape model:
            obj = readObj(obj_file);
            self.faces = obj.f;
            self.verts_b = obj.v;
            if nargin == 7
                self.verts_b = scale*self.verts_b;
            end
            self.verts_i = inert2body'*self.verts_b';
            
            % Read in gravity model:
            load(grav_file,'Cnm','Snm','mu');
            self.Cnm = Cnm;
            self.Snm = Snm;
            self.mu  = mu;
            [~,radii] = normr(self.verts_b);
            self.radius = mean(radii);
            self.radius_max = max(radii);
            self.gravityField = GravityField(self.radius, mu, Cnm, Snm);
            
            % Calculate landmark locations:
            idx = randperm(size(self.verts_b,1));
            self.lmks_b      = self.verts_b(idx(1:num_lmks),:)';
            self.lmk_norms_b = normc(self.lmks_b);
            self.lmks_i      = inert2body'*self.lmks_b;
            self.lmk_norms_i = inert2body'*self.lmk_norms_b;
            
            % Store other data:
            self.rotation_axis = self.inert2body(3,:)';
            self.rotation_rate = rotation_rate;
            self.sun_vec = sun_vec;
            
            % Initialize labels and observation trackers:
            self.lmks_obs = false(1,size(self.lmks_i,2));
            self.lmks_lbl = zeros(1,size(self.lmks_i,2));
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Update parameters about the asteroid:
        function [self] = update(self,dt)
            rotmat = self.rotate(dt);
            self.inert2body = self.inert2body*rotmat;
            self.lmks_i       = self.inert2body'*self.lmks_b;
            self.lmk_norms_i = self.inert2body'*self.lmk_norms_b;
            self.verts_i = self.inert2body'*self.verts_b';
        end
        
        % Get rotation matrix for single time step:
        function [rotmat] = rotate(self,dt,rot_ax,rot_rate)
            if nargin == 2
                rotmat = aa2rotmat(self.rotation_axis, self.rotation_rate*dt);
            else
                rotmat = aa2rotmat(rot_ax, rot_rate*dt);
            end
        end
        
        % Check if observations contain a newly detected lmk:
        function [new_detection,new_inds] = checkForNewDetections(self,visible)
            lmks_obs_prev = self.lmks_obs;
            lmks_obs_new  = self.lmks_obs | visible;
            lmks_obs_sum  = lmks_obs_new + lmks_obs_prev;
            
            % Previously seen lmks will = 2, while new will = 1
            new_inds = lmks_obs_sum == 1;
            if sum(new_inds) > 0
                new_detection = true;
            else
                new_detection = false;
            end
            
            % Store for future:
            self.lmks_obs = lmks_obs_new;
        end
        
        % Create labels for newly detected lmks:
        function [self] = createLabels(self,new_inds,num_tracking)
            current_num = max(self.lmks_lbl);
            self.lmks_lbl(new_inds) = (current_num+1):num_tracking;
        end
    end
    
    %% Public Methods for Visualizations
    methods (Access = public)
        % Draw the asteroid body (useful for animations):
        function [] = drawBody(self)
            if ~self.plotted_body
                self.ptch = patch('Faces',self.faces.v,'Vertices',self.verts_i',...
                                  'FaceColor',1*[1 1 1],'EdgeColor','None',...
                                  'FaceLighting','gouraud','AmbientStrength',0.1,...
                                  'SpecularStrength',0,'FaceAlpha',.5);
                self.ax = gca;
                self.lght = light(self.ax,'Position',1*self.sun_vec);
                hold on; axis equal; rotate3d on
                self.plotted_body = true;
            else
                set(self.ptch,'Vertices',self.verts_i')
            end
        end
        
        % Draw the body in a standalone fashion (useful for a single plot):
        function [p,l] = drawBodyStandalone(self,axs)
            p = patch(axs,'Faces',self.faces.v,'Vertices',self.verts_i',...
                           'FaceColor',1*[1 1 1],'EdgeColor','None',...
                           'FaceLighting','gouraud','AmbientStrength',0.1,...
                           'SpecularStrength',0);
            l = light(axs,'Position',1*self.sun_vec);
            hold on; axis equal; rotate3d on
        end
        
        % Function to draw landmarks:
        function [] = drawLmks(self,visible,varargin)
            if numel(visible) == 1
                inds = true(1,size(self.lmks_i,2));
            elseif length(visible) == size(self.lmks_i,2)
                inds = visible;
            else
                error('Invalid input for VISIBLE')
            end
            if ~self.plotted_lmks
                self.lmk_vis = plot3(self.lmks_i(1,inds),...
                                     self.lmks_i(2,inds),...
                                     self.lmks_i(3,inds),'.g',varargin{:}); hold on
                self.lmk_inv_obs   = plot3(nan,nan,nan,'.b',varargin{:});
                self.lmk_inv_unobs = plot3(self.lmks_i(1,~inds),...
                                           self.lmks_i(2,~inds),...
                                           self.lmks_i(3,~inds),'.r',varargin{:});
                self.plotted_lmks = true;
                axis equal
            else
                set(self.lmk_vis,'XData',self.lmks_i(1,inds),...
                                 'YData',self.lmks_i(2,inds),...
                                 'ZData',self.lmks_i(3,inds));
                inds = ~inds;
                inds_obs   = inds & self.lmks_obs;
                inds_unobs = ~inds_obs & inds;
                set(self.lmk_inv_obs,'XData',self.lmks_i(1,inds_obs),...
                                     'YData',self.lmks_i(2,inds_obs),...
                                     'ZData',self.lmks_i(3,inds_obs));
                set(self.lmk_inv_unobs,'XData',self.lmks_i(1,inds_unobs),...
                                       'YData',self.lmks_i(2,inds_unobs),...
                                       'ZData',self.lmks_i(3,inds_unobs));
            end
        end
        
        % Function to draw estimated (passed in) landmarks:
        function [] = drawLmks_hat(self,X_hat,varargin)
            lmks_h = reshape(X_hat,[],3)';
            if ~self.plotted_lmks_hat
                self.lmks_hat = plot3(lmks_h(1,:),...
                                      lmks_h(2,:),...
                                      lmks_h(3,:),varargin{:}); hold on
                self.plotted_lmks_hat = true;
                axis equal
            else
                set(self.lmks_hat,'XData',lmks_h(1,:),...
                                  'YData',lmks_h(2,:),...
                                  'ZData',lmks_h(3,:));
            end
        end
    end
end