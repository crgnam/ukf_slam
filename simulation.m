clear; matlabrc; clc; close all
addpath(genpath('src'))
addpath(genpath('ukf'))
rng(1)

% Setup:
N = 100; % Number of features to track
dt = 60;
tspan = dt:dt:86400;

f = .055; %(m) Camera focal length
sx = 3; %sensor x in mm
sy = 2; %sensor y in mm
fov = [sx/1000 sy/1000]/f;
K = [f 0 0 0;
     0 f 0 0;
     0 0 1 0];

% Generate fake map:
az = 2*pi*rand(1,N);
el = -pi/2 + pi*rand(1,N);
r = 250*ones(1,N); %(meters) Radius of Bennu
[x,y,z] = sph2cart(az,el,r);
colors = rand(length(x),3); % Colors of each feature (unique)

% Setup orbital parameters:
G = 6.67430*10^-11; %(m3⋅kg–1⋅s–2) Gravitational Constant
M = (7.329*10^10); %(kg) Mass of Bennu
mu = G*M; % Standard gravitational constant of Bennu
a = 700; %(meters) Orbital radius
e = .3;
inc = 45;
w = 0; RAAN = 0; M0 = 0;
[r, v] = orbitalElements2PosVel(a,e,inc,w,RAAN,M0,mu);
apo = a*(1+e); 

% Preallocate memory:
X = zeros(6,length(tspan));
X(:,1) = [r; v];

% Initialize UKF:
r_hat = r;
v_hat = v;

%% Run Filter:
% For the purposes here of this test, it is assumed that ALL features are
% ALWAYS visible
for ii = 1:length(tspan)
    % Simulate the orbit:
    X(:,ii+1) = rk4(@orbitalDynamics,dt,X(:,ii),mu);
    r = X(1:3,ii+1);
    v = X(4:6,ii+1);
    
    % Define a nadir pointing attitude:
    xAxis = -v'/norm(v);
    zAxis = r'/norm(r);
    yAxis = cross(xAxis,zAxis); yAxis = yAxis/norm(yAxis);
    xAxis = cross(yAxis,zAxis); xAxis = xAxis/norm(xAxis);
    rotMat = [xAxis; yAxis; zAxis];
    
    % TEMPORARY: SHOW THE PLOT:
end

%% Animate the Simulation:
figure(1)
subplot(1,2,1)
    scatter3(x,y,z,10,colors,'filled'); hold on; axis equal; grid on
    title('3d simulation')
    xlim([-apo apo])
    ylim([-apo apo])
    zlim([-apo apo])
    ii = 1;
    scale = 100;
    satX = plot3(X(1,ii)+[0 rotMat(1,1)],...
                 X(2,ii)+[0 rotMat(1,2)],...
                 X(3,ii)+[0 rotMat(1,3)],'r');
    satY = plot3(X(1,ii)+[0 rotMat(2,1)],...
                 X(2,ii)+[0 rotMat(2,2)],...
                 X(3,ii)+[0 rotMat(2,3)],'g');
    satZ = plot3(X(1,ii)+[0 rotMat(3,1)],...
                 X(2,ii)+[0 rotMat(3,2)],...
                 X(3,ii)+[0 rotMat(3,3)],'b');
             
subplot(1,2,2)
    image = scatter(x,y,10,colors,'filled'); hold on; grid on; axis equal
    xlim([-fov(1) fov(1)])
    ylim([-fov(2) fov(2)])

% Generate plot:
for ii = 1:length(tspan)
    % Calculate new attitude matrix:
    r = X(1:3,ii);
    v = X(4:6,ii);
    xAxis = -v'/norm(v);
    zAxis = r'/norm(r);
    yAxis = cross(xAxis,zAxis); yAxis = yAxis/norm(yAxis);
    xAxis = cross(yAxis,zAxis); xAxis = xAxis/norm(xAxis);
    rotMat = [xAxis; yAxis; zAxis];
    
    % Update the collected image:
    [imagePoints,inds] = camera(rotMat,r,K,fov, [x;y;z]);
    set(image,'XData',imagePoints(1,:),'YData',imagePoints(2,:),'CData',colors(inds,:))
    
    % Update the pose in the plot:
    rotMat = scale*rotMat;
    set(satX,'XData',X(1,ii)+[0 rotMat(1,1)],...
             'YData',X(2,ii)+[0 rotMat(1,2)],...
             'ZData',X(3,ii)+[0 rotMat(1,3)]);
    set(satY,'XData',X(1,ii)+[0 rotMat(2,1)],...
             'YData',X(2,ii)+[0 rotMat(2,2)],...
             'ZData',X(3,ii)+[0 rotMat(2,3)]);
    set(satZ,'XData',X(1,ii)+[0 rotMat(3,1)],...
             'YData',X(2,ii)+[0 rotMat(3,2)],...
             'ZData',X(3,ii)+[0 rotMat(3,3)]);
    
    drawnow
end