function [handles] = drawOrientation(varargin)
    % h = drawOrientation(h,r,rotMat,scale,___)
    if isstruct(varargin{1})
        handles = varargin{1};
        r = varargin{2};
        rotMat = varargin{3};
        if nargin > 3
            if isnumeric(varargin{4})
                scale = varargin{4};
            end
        else
            scale = 1;
        end
        rotMat = scale*rotMat;
        set(handles.xAxis,'XData',r(1)+[0 rotMat(1,1)],...
                       'YData',r(2)+[0 rotMat(1,2)],...
                       'ZData',r(3)+[0 rotMat(1,3)])
        set(handles.yAxis,'XData',r(1)+[0 rotMat(2,1)],...
                       'YData',r(2)+[0 rotMat(2,2)],...
                       'ZData',r(3)+[0 rotMat(2,3)]);
        set(handles.zAxis,'XData',r(1)+[0 rotMat(3,1)],...
                       'YData',r(2)+[0 rotMat(3,2)],...
                       'ZData',r(3)+[0 rotMat(3,3)]);
    else
        varargin2 = {};
        r = varargin{1};
        rotMat = varargin{2};
        if nargin > 2
            if isnumeric(varargin{3})
                scale = varargin{3};
                if nargin > 3
                    varargin2 = varargin{4:end};
                end
            else
                varargin2 = varargin{3:end};
            end
        else
            scale = 1;
        end
        rotMat = scale*rotMat;
        handles.xAxis = plot3(r(1)+[0 rotMat(1,1)],...
                              r(2)+[0 rotMat(1,2)],...
                              r(3)+[0 rotMat(1,3)],'r',varargin2{:}); hold on
        handles.yAxis = plot3(r(1)+[0 rotMat(2,1)],...
                              r(2)+[0 rotMat(2,2)],...
                              r(3)+[0 rotMat(2,3)],'g',varargin2{:});
        handles.zAxis = plot3(r(1)+[0 rotMat(3,1)],...
                              r(2)+[0 rotMat(3,2)],...
                              r(3)+[0 rotMat(3,3)],'b',varargin2{:});
    end
end