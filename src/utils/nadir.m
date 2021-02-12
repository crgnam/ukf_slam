function [rotMat] = nadir(r,v)
    xAxis = -v'/norm(v);
    zAxis = r'/norm(r);
    yAxis = cross(xAxis,zAxis); yAxis = yAxis/norm(yAxis);
    xAxis = cross(yAxis,zAxis); xAxis = xAxis/norm(xAxis);
    rotMat = [xAxis; yAxis; zAxis];
end