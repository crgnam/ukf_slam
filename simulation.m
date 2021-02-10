clear; matlabrc; clc; close all
addpath(genpath('src'))
addpath(genpath('ukf'))
rng(1)

% Setup:
N = 100; % Number of features to track
dt = 60;
tspan = dt:dt:2*86400;

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
r = 250*ones(1,N) + 20*randn(1,N); %(meters) Radius of Bennu
[x,y,z] = sph2cart(az,el,r);
colors = rand(length(x),3); % Colors of each feature (unique)
% Colors (blue->red based on latitude)
colors(:,1) = 1 - (z+25)/500;
colors(:,2) = 0;
colors(:,3) = (z+250)/500;

% Setup orbital parameters:
G = 6.67430*10^-11; %(m3⋅kg–1⋅s–2) Gravitational Constant
M = (7.329*10^10); %(kg) Mass of Bennu
mu = G*M; % Standard gravitational constant of Bennu
a = 700; %(meters) Orbital radius
e = 0.3;
inc = 45;
w = 0; RAAN = 0; M0 = 0;
[r, v] = orbitalElements2PosVel(a,e,inc,w,RAAN,M0,mu);
apo = a*(1+e); 

% Preallocate memory:
X = zeros(6,length(tspan));
X(:,1) = [r; v];

% Initialize UKF:
P = diag([10*ones(1,3),...
          0.1*ones(1,3),...
          250*ones(1,3*N)]); % Initial State Covariance Matrix
Q = diag([1e-9*ones(1,3),...
          1e-9*ones(1,3),...
          1e-9*ones(1,3*N)]); % Process Noise Covariance Matrix

% TODO: Hardocded measurement covariance matrix, assuming ALL features are
% visible at every time step (valid for first pass)
R = diag(.0001*ones(1,2*N));

% Initial estimate of the map:
rotMat = nadir(r,v);
imagePoints = camera(rotMat,r,K,fov, [x;y;z]);
rays = generateRays(imagePoints,rotMat,K);
projectedPoints = r + norm(r)*rays;

% Initialize estimates:
X_hat = zeros(6 + 3*N, length(tspan));
X_hat(:,1) = [r; v;
              projectedPoints(1,:)';
              projectedPoints(2,:)';
              projectedPoints(3,:)'];
          
% Initialize UKF model arguments:
dynamics_args = {mu};

% UKF tuning parameters:
alpha = 1;
beta  = 10;
kappa = 2;

%% Run Filter:
% For the purposes here of this test, it is assumed that ALL features are
% ALWAYS visible
for ii = 1:length(tspan)
    % Simulate the orbit:
    X(:,ii+1) = rk4(@orbitalDynamics,dt,X(:,ii),mu);
    r = X(1:3,ii+1);
    v = X(4:6,ii+1);
    
    % Define a nadir pointing attitude:
    rotMat = nadir(r,v);
    
    % Calculate the measurement:
    imagePoints = camera(rotMat,r,K,fov, [x;y;z]); %(CURRENTLY NO NOISE)
    measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    
    % Formulate UKF inputs:
    measAvails = true(size(measurement));
    measurement_args = {rotMat,K};
    ukf_args = {dynamics_args, measurement_args};
    
    % Run the UKF:
    [X_hat(:,ii+1), P] = ukf(@ukfOrbit, @ukfCamera,...
                             X_hat(:,ii), dt,...
                             P, Q, R, measAvails, measurement,...
                             alpha, beta, kappa, ukf_args{:});
%     sig3(:,ii+1) = 3*sqrt(diag(P));
end

%% Show the Final Estimated Map:   
estimated_map = reshape(X_hat(7:end,end),[],3)';

figure()
    scatter3(x,y,z,30,'b','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             30,'r','x'); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(X_hat(1,:),X_hat(2,:),X_hat(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(X_hat(1,1),X_hat(2,1),X_hat(3,1),'.r','MarkerSize',20)
    legend('truth map','estimated map','truth trajectory','estimated trajectory')

% % Make animation:
% v = VideoWriter('slam_results.mp4','MPEG-4');
% v.Quality = 90;
% open(v)
% for ii = 1:360
%     view([ii 20])
%     drawnow
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% close(v)

%% Show the history of the SLAM estimates:
% estimated_map = reshape(X_hat(7:end,1),[],3)';
% 
% figure()
% scatter3(x,y,z,30,'b','filled'); hold on; axis equal; grid on
% est_map = scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
%          30,'r','x'); hold on; axis equal; grid on
% true_sat = plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20);
% est_sat = plot3(X_hat(1,1),X_hat(2,1),X_hat(3,1),'.r','MarkerSize',20);
% legend('truth map','estimated map','truth trajectory','estimated trajectory')
%     
%     % Make animation:
%     v = VideoWriter('slam_results.mp4','MPEG-4');
%     v.Quality = 90;
%     open(v)
%     for ii = 1:size(X_hat,2)
%         estimated_map = reshape(X_hat(7:end,ii),[],3)';
%         set(true_sat,'XData',X(1,ii),'YData',X(2,ii),'ZData',X(3,ii));
%         set(est_sat,'XData',X_hat(1,ii),'YData',X_hat(2,ii),'ZData',X_hat(3,ii));
%         set(est_map,'XDAta',estimated_map(1,:),'YData',estimated_map(2,:),'ZData',estimated_map(3,:))
%         view([ii*(360/size(X_hat,2)) 20])
%         drawnow
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     end
%     close(v)

%% Animate the Simulation:
% figure(1)
% subplot(1,2,1)
%     scatter3(x,y,z,30,colors,'filled'); hold on; axis equal; grid on
%     title('3d simulation')
%     xlim([-apo apo])
%     ylim([-apo apo])
%     zlim([-apo apo])
%     ii = 1;
%     scale = 100;
%     satX = plot3(X(1,ii)+[0 rotMat(1,1)],...
%                  X(2,ii)+[0 rotMat(1,2)],...
%                  X(3,ii)+[0 rotMat(1,3)],'r');
%     satY = plot3(X(1,ii)+[0 rotMat(2,1)],...
%                  X(2,ii)+[0 rotMat(2,2)],...
%                  X(3,ii)+[0 rotMat(2,3)],'g');
%     satZ = plot3(X(1,ii)+[0 rotMat(3,1)],...
%                  X(2,ii)+[0 rotMat(3,2)],...
%                  X(3,ii)+[0 rotMat(3,3)],'b');
%     ax = gca;
%     ax.Projection = 'perspective';
%              
% subplot(1,2,2)
%     image = scatter(x,y,30,colors,'filled'); hold on; grid on; axis equal
%     xlim([-fov(1) fov(1)])
%     ylim([-fov(2) fov(2)])
%     plot([0 fov(1)/5],[0 0],'r')
%     plot([0 0],[0 fov(1)/5],'g')
% 
% % Generate plot:
% for ii = 1:length(tspan)
%     % Calculate new attitude matrix:
%     r = X(1:3,ii);
%     v = X(4:6,ii);
%     rotMat = nadir(r,v);
%     
%     % Update the collected image:
%     [imagePoints,inds] = camera(rotMat,r,K,fov, [x;y;z]);
%     set(image,'XData',imagePoints(1,:),'YData',imagePoints(2,:),'CData',colors(inds,:))
%     
%     % Update the pose in the plot:
%     rotMat = scale*rotMat;
%     set(satX,'XData',X(1,ii)+[0 rotMat(1,1)],...
%              'YData',X(2,ii)+[0 rotMat(1,2)],...
%              'ZData',X(3,ii)+[0 rotMat(1,3)]);
%     set(satY,'XData',X(1,ii)+[0 rotMat(2,1)],...
%              'YData',X(2,ii)+[0 rotMat(2,2)],...
%              'ZData',X(3,ii)+[0 rotMat(2,3)]);
%     set(satZ,'XData',X(1,ii)+[0 rotMat(3,1)],...
%              'YData',X(2,ii)+[0 rotMat(3,2)],...
%              'ZData',X(3,ii)+[0 rotMat(3,3)]);
%     drawnow
% end