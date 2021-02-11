clear; matlabrc; clc; close all; rng(1);
import = @(x) addpath(genpath(x));
import('data')
import('lib')
import('plotting')
import('src')
import('ukf')

%% Setup:
% Simulation Setup:
N = 50; % number of features to be simulated and tracked
dt = 60;
tspan = dt:dt:1*86400;


% Setup camera (pinhole model):
f = .055; %(m) Camera focal length
res_x = 1920; %(pix)
res_y = 1080; %(pix)
sx = 3/1000; %(m) sensor x 
sy = 2/1000; %(m) sensor y
fov = [sx sy]/f;
K = [f 0 0 0;
     0 f 0 0;
     0 0 1 0];
 
% Orbit setup:
semi_major   = 700; %(meters) Orbital radius
eccentricity = 0.3;
inclination  = 45;

%% Initialize Simulation:
% Generate simulated landmarks:
bennu = Bennu(N);

% Calculate orbital state vector:
[r,v] = orbitalElements2PosVel(semi_major,eccentricity,inclination,0,0,0,bennu.mu);
apoapsis = 2*semi_major*(1+eccentricity); 
tspan = dt:dt:2*pi*sqrt(semi_major^3/bennu.mu);

% Preallocate memory:
X = zeros(6,length(tspan));

% Assign initial values:
X(:,1) = [r; v];

%% Initialize Filters:
% UKF tuning parameters:
alpha = 1;
beta  = 10;
kappa = 2;

% Initial State Covariance Matrix:
P = diag([10*ones(1,3),...
          0.1*ones(1,3),...
          250*ones(1,3*N)]);
      
% Process Noise Covariance Matrix:
Q = diag([1e-9*ones(1,3),...
          1e-9*ones(1,3),...
          1e-9*ones(1,3*N)]);

% Measurement Covariance Matrix:
sig_meas = fov(1)/res_x;
R = diag((sig_meas^2)*ones(1,2*N));

% Initial estimate of the map:
rotMat = nadir(r,v);
imagePoints = camera(rotMat,r,K,fov, bennu.lmks,nan, sig_meas);
rays = generateRays(imagePoints,rotMat,K);
projectedPoints = r + norm(r)*rays;

% Initialize estimates:
X_hat = zeros(6 + 3*N, length(tspan));
X_hat(:,1) = [r; v;
              projectedPoints(1,:)';
              projectedPoints(2,:)';
              projectedPoints(3,:)'];
          
% Initialize UKF model arguments:
dynamics_args = {bennu.mu};

% 3-sigma bound init:
sig3 = zeros(size(P,1),length(tspan));
sig3_rts = sig3;
sig3(:,1) = 3*sqrt(diag(P));
avails_history = zeros(2*N,length(tspan));
meas_history   = zeros(2*N,length(tspan));
rotMat_hist    = zeros(3,3,length(tspan));
rotMat_hist(:,:,1) = rotMat;

% Smoother initialization:
P_hist = zeros(size(P,1),size(P,2),length(tspan));
P_hist(:,:,1) = P;
disp('Completed Initialization')

%% Run Simulation:
visualize  = false;
save_video = false;
if visualize
    figure('units','normalized','outerposition',[0 0 1 1])
    bennu.drawBody()
    camva(2)
    cam = Attitude(150,'LineWidth',2);
    traj = plot3(r(1),r(2),r(3),'r');
    grid on
    view([45 20])
    xlim([-apoapsis apoapsis])
    ylim([-apoapsis apoapsis])
    zlim([-apoapsis apoapsis])
end
if save_video
    vid = VideoWriter('animation.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end

for ii = 1:length(tspan)
    % Simulate the orbit:
    X(:,ii+1) = rk4(@orbitalDynamics,dt,X(:,ii),bennu.mu);
    r = X(1:3,ii+1);
    v = X(4:6,ii+1);
    
    % Define a nadir pointing attitude:
    rotMat = nadir(r,v);
    
    % Calculate the measurement:
    [imagePoints,visible] = camera(rotMat,r,K,fov, bennu.lmks,bennu.lmk_norms, sig_meas);
    measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    
    % Visualize:
    if visualize 
        bennu.drawLmks(visible,'MarkerSize',10)
        view([20+ii*360/length(tspan),20])
        cam.draw(r,rotMat)
        set(traj,'XData',X(1,1:ii),'YData',X(2,1:ii),'ZData',X(3,1:ii))
        drawnow
    end
    
    % Make animation:
    if save_video
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    
    % Formulate UKF inputs:
    measAvails = [visible visible]';
    measurement_args = {rotMat,K};
    ukf_args = {dynamics_args, measurement_args};
    
    % Run the UKF:
    [X_hat(:,ii+1), P] = ukf(@ukfOrbit, @ukfCamera,...
                             X_hat(:,ii), dt,...
                             P, Q, R, measAvails, measurement,...
                             alpha, beta, kappa, ukf_args{:});
    sig3(:,ii+1) = 3*sqrt(diag(P));
    P_hist(:,:,ii+1) = P;
    avails_history(:,ii+1) = measAvails;
    meas_history(1:sum(measAvails),ii+1) = measurement;
    rotMat_hist(:,:,ii+1) = rotMat;
end
if save_video
    close(vid)
end
disp('Completed Forward Simulation')

%% Show Horns Method 
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Horn's method to transform to truth space:
[regParams,~, errorStats_ukf] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);

save_video = true;
if save_video
    vid = VideoWriter('horns.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'g','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             100,'m','o','LineWidth',2); hold on; axis equal; grid on
    legend('True Map','Estimated Map','location','southeast')
    camva(1)
    xlim([-apoapsis apoapsis])
    ylim([-apoapsis apoapsis])
    zlim([-apoapsis apoapsis])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
if save_video
    for ii = 1:360
        view([ii 20])
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    close(vid)
end

%% Show Kabsch Method:
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Kabsch algorithm to transform to truth space:
[R_ukf, t_ukf] = kabsch(bennu.lmks, estimated_map);
estimated_map = R_ukf'*(estimated_map - t_ukf);

save_video = true;
if save_video
    vid = VideoWriter('kabsch.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'g','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             100,'m','o','LineWidth',2); hold on; axis equal; grid on
    legend('True Map','Estimated Map','location','southeast')
    camva(1)
    xlim([-apoapsis apoapsis])
    ylim([-apoapsis apoapsis])
    zlim([-apoapsis apoapsis])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
if save_video
    for ii = 1:360
        view([ii 20])
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    close(vid)
end