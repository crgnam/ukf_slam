classdef Asteroid < handle
    properties 
        % Shape model:
        faces
        verts
        
        % Landmarks:
        lmks
        lmk_norms
        
        % Environment:
        rotmat
        radius
        rotation_axis
        rotation_rate
        sun_vec
        
        % Gravity parameters:
        gravityField
        Cnm
        Snm
        mu
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = Asteroid(obj_file, grav_file, orientation,...
                                   rotation_rate,sun_vec, num_lmks, scale)
            % Get orientation:
            self.rotmat = orientation;
            
            % Read in shape model:
            obj = readObj(obj_file);
            self.faces = obj.f;
            self.verts = obj.v;
            if nargin == 7
                self.verts = scale*self.verts;
            end
            
            % Read in gravity model:
            load(grav_file,'Cnm','Snm','mu');
            self.Cnm = Cnm;
            self.Snm = Snm;
            self.mu  = mu;
            [~,radii] = normr(self.verts);
            self.radius = mean(radii);
            self.gravityField = GravityField(self.radius, mu, Cnm, Snm);
            
            % Calculate landmark locations:
            idx = randperm(size(self.verts,1));
            self.lmks = self.verts(idx(1:num_lmks),:)';
            self.lmk_norms = normc(self.lmks);
            
            % Store other data:
            self.rotation_axis = self.rotmat(3,:)';
            self.rotation_rate = rotation_rate;
            self.sun_vec = sun_vec;
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Update parameters about the asteroid:
        function [self] = update(self,dt)
            rotmat2 = aa2rotmat(self.rotation_axis, self.rotation_rate*dt);
            self.rotmat = self.rotmat*rotmat2;
        end
    end
end