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
apoapsis = semi_major*(1+eccentricity); 

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
sig3_post = sig3;
sig3(:,1) = 3*sqrt(diag(P));
avails_history = zeros(2*N,length(tspan));
meas_history   = zeros(2*N,length(tspan));
rotMat_hist    = zeros(3,3,length(tspan));
rotMat_hist(:,:,1) = rotMat;

% Smoother initialization:
P_hist = zeros(size(P,1),size(P,2),length(tspan));
P_hist(:,:,1) = P;

%% Run Simulation:
visualize  = false;
save_video = false;
if visualize
    figure('units','normalized','outerposition',[0 0 1 1])
    bennu.drawBody()
    cam = Attitude(150,'LineWidth',2);
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
        cam.draw(r,rotMat)
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

%% Show Forward Results:  
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Kabsch algorithm to transform to truth space:
[regParams,~,~] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R2 = regParams.R;
t2 = regParams.t;
S2 = regParams.s;
estimated_map = (1/S2)*R2'*(estimated_map - t2);
M_ukf = (1/S2)*R2'*(X_hat(1:3,:) - t2);

if save_video
    vid = VideoWriter('ukf_results.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'b','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             40,'r','x'); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(M_ukf(1,:),M_ukf(2,:),M_ukf(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(M_ukf(1,1),M_ukf(2,1),M_ukf(3,1),'.r','MarkerSize',20)
    legend('True Map','Estimated Map','True Trajectory','Estimated Trajectory',...
           'location','south')
    title('Forward UKF Results')
if save_video
    for ii = 1:360
        view([ii 20])
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    close(vid)
end
    
%% Run the smoother:
[M,P_hist,D] = urts(X_hat,P_hist,@ukfOrbit,dt,Q,dynamics_args,alpha,beta,kappa,0,1);

% Get the 3-sigma bounds:
% for ii = 1:size(P_hist,3)
%     sig3_post(:,ii) = 3*sqrt(diag(P_hist(:,:,ii)));
% end
    
%% Show Smoother Results:
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Horn's Method to transform to truth space:
[regParams,Bfit,ErrorStats] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R2 = regParams.R;
t2 = regParams.t;
S2 = regParams.s;

% Transform map and states:
estimated_map = (1/S2)*R2'*(estimated_map - t2);
M_rts(1:3,:) = (1/S2)*R2'*(M(1:3,:) - t2);

if save_video
    vid = VideoWriter('urts_results.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'b','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             40,'r','x'); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(M_rts(1,:),M_rts(2,:),M_rts(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(M_rts(1,1),M_rts(2,1),M_rts(3,1),'.r','MarkerSize',20)
    legend('True Map','Estimated Map','True Trajectory','Estimated Trajectory',...
           'location','south')
    title('Unscented RTS Smoother Results')
    if save_video
        for ii = 1:360
            view([ii 20])
            frame = getframe(gcf);
            writeVideo(vid,frame);
        end
        close(vid)
    end
        
%% Calculate Residuals:
pt_resid = zeros(size(meas_history));
pf_resid = zeros(size(meas_history));
for ii = 1:size(X_hat,2)
    % Get the measurement availability:
    avails = avails_history(:,ii)==1;
    
    % Calculate pass-through predicted measurement:
    pt_predict = ukfCamera(dt, X_hat(:,ii),rotMat_hist(:,:,ii),K, avails);
    pf_predict = ukfCamera(dt, M(:,ii),    rotMat_hist(:,:,ii),K, avails);
    
    % Calculate the measurement residuals:
    pt_resid(1:sum(avails),ii) = pt_predict - meas_history(1:sum(avails),ii);
    pf_resid(1:sum(avails),ii) = pf_predict - meas_history(1:sum(avails),ii);
end

% Plot the residuals:
pt_pixel = pt_resid(1:N,:);
pt_row   = pt_resid(N+1:end,:);
pf_pixel = pf_resid(1:N,:);
pf_row   = pf_resid(N+1:end,:);

figure()
    subplot(1,2,1)
        plot(pt_pixel*res_x,'.b'); hold on
        plot(pt_row*res_y,'.r');
        ylim([-5 5])
        grid on

    subplot(1,2,2)
        plot(pf_pixel*res_x,'.b'); hold on
        plot(pf_row*res_y,'.r');
        ylim([-5 5])
        grid on