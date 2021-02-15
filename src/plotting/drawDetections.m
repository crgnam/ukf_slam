function [] = drawDetections(t_plt,detection_hist,y,varargin)
    t_plt = t_plt(1,1:end);
    for ii = 1:length(detection_hist)
        t = t_plt(detection_hist(ii));
        plot([t t],[y(1) y(2)],varargin{:}); hold on
    end
end