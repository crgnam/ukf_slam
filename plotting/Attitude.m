classdef Attitude < handle
    properties (Access = public)
        scale
        xAxis
        yAxis
        zAxis
    end
    
    %% Constructor
    methods (Access = public)
        function [self] = Attitude(scale,varargin)
            r = nan(3,1);
            rotMat = nan(3);
            self.scale = scale;
            self.xAxis = plot3(r(1)+[0 rotMat(1,1)],...
                               r(2)+[0 rotMat(1,2)],...
                               r(3)+[0 rotMat(1,3)],'r',varargin{:});
            self.yAxis = plot3(r(1)+[0 rotMat(2,1)],...
                               r(2)+[0 rotMat(2,2)],...
                               r(3)+[0 rotMat(2,3)],'g',varargin{:});
            self.zAxis = plot3(r(1)+[0 rotMat(3,1)],...
                               r(2)+[0 rotMat(3,2)],...
                               r(3)+[0 rotMat(3,3)],'b',varargin{:});
        end
    end
    
    %% Public Methods:
    methods (Access = public)
        function [self] = draw(self,r,rotMat)
            rotMat = self.scale*rotMat;
            set(self.xAxis,'XData',r(1)+[0 rotMat(1,1)],...
                           'YData',r(2)+[0 rotMat(1,2)],...
                           'ZData',r(3)+[0 rotMat(1,3)]);
            set(self.yAxis,'XData',r(1)+[0 rotMat(2,1)],...
                           'YData',r(2)+[0 rotMat(2,2)],...
                           'ZData',r(3)+[0 rotMat(2,3)]);
            set(self.zAxis,'XData',r(1)+[0 rotMat(3,1)],...
                           'YData',r(2)+[0 rotMat(3,2)],...
                           'ZData',r(3)+[0 rotMat(3,3)]);
        end
    end
end