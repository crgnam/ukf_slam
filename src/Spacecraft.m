classdef Spacecraft < handle
    properties
        % States of the spacecraft:
        r
        v
        rotmat
        
        % camera object:
        camera
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
            X = rk4(@self.dynamics,dt,[self.r; self.v],body);
            self.r = X(1:3);
            self.v = X(4:6);
        end
        
        % Take a measurement given a body with a set of lmks:
        function [lmks,visible] = image(self,body)
            
        end
    end
    
    methods (Access = private)
        function [dX] = dynamics(~,~,X,body)
            a = body.gravityField.acceleration(X(1:3), body.rotmat);
            dX = [X(4:6); a];
        end
    end
end