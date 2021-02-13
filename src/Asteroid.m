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
        
        % Environment:
        inert2body
        radius
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
        lmk_inv
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
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Update parameters about the asteroid:
        function [self] = update(self,dt)
            rotmat = aa2rotmat(self.rotation_axis, self.rotation_rate*dt);
            self.inert2body = self.inert2body*rotmat;
            self.lmks_i       = self.inert2body'*self.lmks_b;
            self.lmk_norms_i = self.inert2body'*self.lmk_norms_b;
            self.verts_i = self.inert2body'*self.verts_b';
        end
        
        function [] = drawBody(self)
            if ~self.plotted_body
                self.ptch = patch('Faces',self.faces.v,'Vertices',self.verts_i',...
                                  'FaceColor',1*[1 1 1],'EdgeColor','None',...
                                  'FaceLighting','gouraud','AmbientStrength',0.1,...
                                  'SpecularStrength',0);
                self.ax = gca;
                self.lght = light(self.ax,'Position',1*self.sun_vec);
                hold on; axis equal; rotate3d on
                self.plotted_body = true;
            else
                set(self.ptch,'Vertices',self.verts_i')
            end
%             set(gca,'Color','k')
        end
        
        % Draw the body in a standalone fashion:
        function [p,l] = drawBodyStandalone(self,axs)
            p = patch(axs,'Faces',self.faces.v,'Vertices',self.verts,...
                           'FaceColor',[0.5 0.5 0.5],'EdgeColor','None',...
                           'FaceLighting','gouraud','AmbientStrength',0.5,...
                           'SpecularStrength',0);
            l = light(axs,'Position',1*self.sun_vec);
            hold on; axis equal
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
                                     self.lmks_i(3,inds),'.g',varargin{:});
                self.lmk_inv = plot3(self.lmks_i(1,~inds),...
                                     self.lmks_i(2,~inds),...
                                     self.lmks_i(3,~inds),'.r',varargin{:});
                self.plotted_lmks = true;
                axis equal
            else
                set(self.lmk_vis,'XData',self.lmks_i(1,inds),...
                                 'YData',self.lmks_i(2,inds),...
                                 'ZData',self.lmks_i(3,inds));
                set(self.lmk_inv,'XData',self.lmks_i(1,~inds),...
                                 'YData',self.lmks_i(2,~inds),...
                                 'ZData',self.lmks_i(3,~inds));
            end
        end
    end
end