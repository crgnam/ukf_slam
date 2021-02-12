classdef Asteroid < handle
    properties 
        % Shape model:
        faces
        verts
        
        % Landmarks:
        lmks
        lmk_norms
        
        % Environment:
        sun_vec
        
        % Gravity parameters:
        Cnm
        Snm
        mu
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = Asteroid(obj_file, grav_file, sun_vec, num_lmks, scale)
            % Read in shape model:
            obj = readObj(obj_file);
            self.faces = obj.f;
            self.verts = obj.v;
            if nargin == 5
                self.verts = scale*self.verts;
            end
            
            % Read in gravity model:
            load(grav_file,'Cnm','Snm','mu');
            self.Cnm = Cnm;
            self.Snm = Snm;
            self.mu  = mu;
            
            % Calculate landmark locations:
            idx = randperm(size(self.verts,1));
            self.lmks = self.verts(idx(1:num_lmks),:)';
            self.lmk_norms = normc(self.lmks);
            
            % Store other data:
            self.sun_vec = sun_vec;
        end
    end
    
    %% Public Methods
    methods (Access = public)
        
    end
end