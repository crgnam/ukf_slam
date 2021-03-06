function [handles] = drawOrientation(r,rotMat,lineType,varargin)
    handles.xAxis = plot3(r(1)+[0 rotMat(1,1)],...
                          r(2)+[0 rotMat(1,2)],...
                          r(3)+[0 rotMat(1,3)],'r','LineStyle',lineType,varargin{:}); hold on
    handles.yAxis = plot3(r(1)+[0 rotMat(2,1)],...
                          r(2)+[0 rotMat(2,2)],...
                          r(3)+[0 rotMat(2,3)],'g','LineStyle',lineType,varargin{:});
    handles.zAxis = plot3(r(1)+[0 rotMat(3,1)],...
                          r(2)+[0 rotMat(3,2)],...
                          r(3)+[0 rotMat(3,3)],'b','LineStyle',lineType,varargin{:});
end