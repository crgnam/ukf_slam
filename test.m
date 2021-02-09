clear; matlabrc; clc; close all;
addpath(genpath('src'))

worldPoints = randn(3,100);
colors = rand(length(worldPoints),3);

rotMat = e2a(0,45,180);
position = [0; 0; 0];
% xAxis = [1 0 0];
% zAxis = position'/norm(position);
% yAxis = cross(xAxis,zAxis); yAxis = yAxis/norm(yAxis);
% xAxis = cross(yAxis,zAxis); xAxis = xAxis/norm(xAxis);
% rotMat = [xAxis; yAxis; zAxis];


f = .5;
fov = [3 2];
K = [f 0 0;
     0 f 0;
     0  0 1];

[imagePoints,inds] = camera(rotMat,position,K,fov, worldPoints);

subplot(1,2,1)
    scatter(imagePoints(1,:),imagePoints(2,:),10,colors(inds,:),'filled'); grid on; axis equal
    xlim([-fov(1) fov(1)])
    ylim([-fov(2) fov(2)])

subplot(1,2,2)
    scatter3(worldPoints(1,:)',worldPoints(2,:)',worldPoints(3,:)',10,colors,'filled'); hold on; rotate3d on; axis equal; grid on
    plot3(position(1)+[0 rotMat(1,1)],...
          position(2)+[0 rotMat(1,2)],...
          position(3)+[0 rotMat(1,3)],'r')
    plot3(position(1)+[0 rotMat(2,1)],...
          position(2)+[0 rotMat(2,2)],...
          position(3)+[0 rotMat(2,3)],'g')
    plot3(position(1)+[0 rotMat(3,1)],...
          position(2)+[0 rotMat(3,2)],...
          position(3)+[0 rotMat(3,3)],'b')