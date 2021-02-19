classdef UWB < handle
    properties
        r % 3xN positions of N transceivers
        v % 3xN velocities of N transceivers
        num_transceivers
        
        % Measurement parameters:
        range_unc
        angle_unc
        
        % Animation stuff:
        plotted_positions = false;
        positions_handle
        plotted_traj = false;
        traj_handle
        plotted_connections = false;
        connects_handle
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = UWB(r,v,range_unc,angle_unc)
            self.r = r;
            self.v = v;
            self.num_transceivers = size(self.r,2);
            
            self.range_unc = range_unc;
            self.angle_unc = angle_unc;
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Propagate all of the particles:
        function [self] = propagate(self,dt,body)
            for ii = 1:self.num_transceivers
                X = rk4(@gravFieldDynamics,dt,[self.r(:,ii); self.v(:,ii)],body);
                self.r(:,ii) = X(1:3);
                self.v(:,ii) = X(4:6);
            end
        end
        
        % Calculate all of the range measurements:
        function [measurement,avail,vec_pairs] = measureRanges(self,satellite,body,std_meas)
            % Generate all pairs of connections:
            vec_pairs = combvec(self.r);
            vec_pairs = [vec_pairs;
                         repmat(satellite.r',self.num_transceivers,1), self.r';
                         self.r', repmat(satellite.r',self.num_transceivers,1)];

            % Find pairs which are blocked by the body:
            origins = vec_pairs(:,1:3);
            [rays,ranges] = normr(vec_pairs(:,4:6) - origins);
            lines  = [origins rays];
            sphere = [0 0 0 body.radius];
            intersects = intersectLineSphere(lines, sphere);
            intersects = sum([intersects(1:size(intersects,1)/2,:),...
                              intersects(1+size(intersects,1)/2:end,:)],2);
                 
            % Create boolean array of available measurements, and measurements:
            avail = isnan(intersects);
            measurement = ranges(avail) + std_meas*randn(size(ranges(avail),1),1);
        end
        
        % Calculate all of the PDoA angle measurements:
        function [measurement] = measureAngles(self,std_meas)
            
        end
    end
    
    %% Public Methods for Visualization:
    methods (Access = public)
        function [self] = draw(self,varargin)
            if ~self.plotted_positions
                self.positions_handle = plot3(self.r(1,:),self.r(2,:),self.r(3,:),...
                                              varargin{:}); hold on
                self.plotted_positions = true;
            else
                set(self.positions_handle,'XData',self.r(1,:),...
                                          'YData',self.r(2,:),...
                                          'ZData',self.r(3,:))
            end
        end
        
        function [self] = drawTraj(self,r_hist,varargin)
            if ~self.plotted_traj
                for ii = 1:self.num_transceivers
                    self.traj_handle(ii) = plot3(r_hist(1+(ii-1)*3,:),...
                                                 r_hist(2+(ii-1)*3,:),....
                                                 r_hist(3+(ii-1)*3,:),...
                                                 varargin{:}); hold on
                    self.plotted_traj = true;
                end
            else
                for ii = 1:self.num_transceivers
                    set(self.traj_handle(ii),'XDAta',r_hist(1+(ii-1)*3,:),...
                                             'YData',r_hist(2+(ii-1)*3,:),...
                                             'ZData',r_hist(3+(ii-1)*3,:))
                end
            end
        end
        
        function [self] = drawConnections(self,satellite,body,on,varargin)
            % Generate all pairs:
            if on
                [~,avail,vec_pairs] = self.measureRanges(satellite,body,0);
                vec_pairs(~avail,:) = nan;
            else
                num_pairs = factorial(self.num_transceivers)/factorial(self.num_transceivers - 2);
                vec_pairs =  nan(num_pairs+2*self.num_transceivers,6);
            end
            
            % Update the drawing:
            if ~self.plotted_connections
                for ii = 1:size(vec_pairs,1)
                   self.connects_handle(ii) = plot3([vec_pairs(ii,1) vec_pairs(ii,4)],...
                                                    [vec_pairs(ii,2) vec_pairs(ii,5)],...
                                                    [vec_pairs(ii,3) vec_pairs(ii,6)],varargin{:});
                end
                self.plotted_connections = true;
            else
                for ii = 1:size(vec_pairs,1)
                   set(self.connects_handle(ii),'XData',[vec_pairs(ii,1) vec_pairs(ii,4)],...
                                                'YData',[vec_pairs(ii,2) vec_pairs(ii,5)],...
                                                'ZData',[vec_pairs(ii,3) vec_pairs(ii,6)]);
                end
            end
        end
    end
end