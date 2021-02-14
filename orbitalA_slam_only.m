clear; matlabrc; clc; close all; rng(1);
import = @(x) addpath(genpath(x));
import('data')
import('lib')
import('src')
import('filters')

%% Simulation Setup:
% Objects to be tracked:
num_lmks = 50; % Number of features to be tracked

% Orbit setup: (currently using rough orbital A)
semi_major   = 1000; %(meters) Orbital radius
eccentricity = 0.2;
inclination  = 80;
arg_peri     = 0;
right_asc    = 0;
mean_anom    = 0;
apoapsis = semi_major*(1+eccentricity); 

% Bennu states:
sun_vec             = [0;1;0]; % Sun vector (for illumination visibility of lmks)
bennu_inertial2body = ea2rotmat(0,0,0,'321'); % Orientation of Bennu
bennu_w             = 0*2*pi/(4.3*3600); %(rad/s) Rotation rate of bennu

% Setup camera (pinhole model):
f = .055; %(m) Camera focal length
resolution = [1920 1080]; %(pixels) Image dimensions
sensor     = [3 2]/1000; %(m) Sensor dimensions
fov = [sensor(1) sensor(2)]/f;
K = [f 0 0 0; 0 f 0 0; 0 0 1 0]; % Camera projection matrix:
std_meas = fov(1)/resolution(1); % Measurement uncertainty of a single pixel
max_view_angle = deg2rad(70); % Maximum angle away from the normal at which a feature is "visible"
 
% Time setup:
dt = 1*60;
duration = 1*86400;

%% Filter Initialization:
alpha = 1e-4;
beta  = 2000;
kappa = 4;

% Estimation covariance initial values:
p_r_unc = 10;
p_v_unc = 10;
p_lmk_unc = [200; 200; 200]; % Uncertainties in lmks x,y,z

% Process noise covariance initial values:
q_r = 1e-7;
q_v = 1e-7;
q_lmk = [1e-7; 1e-7; 1e-7];


%% Initialize:
% Instantiate camera and asteroid objects:
camera = Camera(K,fov,sensor,resolution,max_view_angle);
bennu  = Asteroid('data/bennu.obj','data/bennugrav.mat',...
                  bennu_inertial2body,bennu_w,sun_vec,num_lmks,1000);

% Calculate important initial values:
tspan = dt:dt:duration;
[r,v] = orb2rv(semi_major,eccentricity,inclination,...
               arg_peri,right_asc,mean_anom,bennu.mu);
                           
% Instantiate the orex vehicle:
rotmat = nadir(r,v);
orex   = Spacecraft(r,v,rotmat,camera);


%% Pre-allocate memory for logging data:
L = length(tspan);

r_hist = zeros(3,L); r_hist(:,1) = r; 
v_hist = zeros(3,L); v_hist(:,1) = v;

lmks_hist = zeros(3,size(bennu.lmks_i,2),L);
visible_hist = zeros(L,num_lmks);

meas_hist   = zeros(2*num_lmks,length(tspan));
avails_hist = zeros(2*num_lmks,length(tspan));

%% Initialize Map Estimate:
r_hat = r + p_r_unc*randn(3,1);
v_hat = v + p_v_unc*randn(3,1);

% Take first image:
[image_lmks,visible] = orex.image(bennu,0);

% Create labels for new features:
[new_detection,new_inds] = bennu.checkForNewDetections(visible);
num_tracking = sum(bennu.lmks_obs);
bennu.createLabels(new_inds,num_tracking);

% Get initial map:
radius_hat      = camera.radius_estimate(r_hat,orex.r,bennu.radius);
projectedPoints = initializeMap(camera,image_lmks,orex.rotmat,r_hat,radius_hat);

% Preallocate estimate:
X_hat = zeros(6+3*num_tracking,L);
sig3 = X_hat;

% Store as initial estimate:
X_hat(:,1) = [r; v; projectedPoints];

% Initial Covariance matrices:
P = diag([p_r_unc*ones(1,3),...
          p_v_unc*ones(1,3),...
          p_lmk_unc(1)*ones(1,num_tracking),...
          p_lmk_unc(2)*ones(1,num_tracking),...
          p_lmk_unc(3)*ones(1,num_tracking)]);
      
Q = diag([q_r*ones(1,3),...
          q_v*ones(1,3),...
          q_lmk(1)*ones(1,num_tracking),...
          q_lmk(2)*ones(1,num_tracking),...
          q_lmk(3)*ones(1,num_tracking)]);
      
R = diag((std_meas^2)*ones(1,2*num_tracking));

sig3(:,1) = 3*sqrt(diag(P));

% Format arguments for ukf input:
dynamics_args    = {bennu};
measurement_args = {orex};
ukf_args = {dynamics_args, measurement_args};

%% See the initial map estimate:
% figure('units','normalized','outerposition',[0 0 1 1])
% orex.drawRays(image_lmks,radius_hat,r_hat);
% bennu.drawBody();
% % bennu.lght = [];
% set(bennu.ptch,'AmbientStrength',.8,'FaceAlpha',0.5)
% bennu.drawLmks(visible,'MarkerSize',30);
% bennu.drawLmks_hat(X_hat(7:end,1),'ob','MarkerSize',10);
% setLims(apoapsis)
% camva(2)

%% Run Simulation:
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    
    % Collect measurement:
    [image_lmks,visible,lmk_inds] = orex.image(bennu,std_meas);
    
    % Generate labels if needed:
    [new_detection,new_inds] = bennu.checkForNewDetections(visible);
    
    % Get initial estimate:
    if new_detection
        % Update the number of all lmks being tracked:
        num_tracking = sum(bennu.lmks_obs);
        
        % Create labels for the newly created objects:
        bennu.createLabels(new_inds,num_tracking);
        
        % Determine which image points correspond to a new feature:
        lmk_new_inds = 1:size(bennu.lmks_i,2);
        lmk_new_inds = lmk_new_inds(new_inds);
        [~,idx] = intersect(lmk_inds,lmk_new_inds);
        
        % Get initial estimate for newly detected feature:
        radius_hat = camera.radius_estimate(X_hat(1:3,ii),orex.r,bennu.radius);
        X_hat_new  = initializeMap(camera,image_lmks(:,idx),...
                                   orex.rotmat,X_hat(1:3,ii),radius_hat);
        
        % Augment filter if needed:
        [X_hat,P,Q,sig3] = ukf_augment(X_hat,P,Q,sig3, ii, X_hat_new,...
                                       p_lmk_unc, q_lmk);
        
        P = diag([p_r_unc*ones(1,3),...
                  p_v_unc*ones(1,3),...
                  p_lmk_unc(1)*ones(1,num_tracking),...
                  p_lmk_unc(2)*ones(1,num_tracking),...
                  p_lmk_unc(3)*ones(1,num_tracking)]);

        Q = diag([q_r*ones(1,3),...
                  q_v*ones(1,3),...
                  q_lmk(1)*ones(1,num_tracking),...
                  q_lmk(2)*ones(1,num_tracking),...
                  q_lmk(3)*ones(1,num_tracking)]);
        % Increase measurement covariance matrix for all possible features
        % being tracked now:
        R = diag((std_meas^2)*ones(1,2*num_tracking));
        fprintf('Added %i new features to the filter at time step %i\n',sum(new_inds),ii)
    end
    
    % Sort measurements by the corresponding landmark label:
    labeled_image_lmks = sortrows([bennu.lmks_lbl(visible); image_lmks]')';
    measurement = [labeled_image_lmks(2,:)'; labeled_image_lmks(3,:)'];
    
    % Determine which measurements of the states are available:
    [~,ia] = intersect(1:num_tracking,bennu.lmks_lbl(visible));
    measAvailBools = false(size(1:num_tracking));
    measAvailBools(ia) = true;
    measAvails = [measAvailBools'; measAvailBools'];
    
    % Run the filter:======================================================
    [X_hat(:,ii+1),P] = ukf(@ukfOrbit, @ukfCamera,...
                            X_hat(:,ii), dt,...
                            P, Q, R, measAvails, measurement,...
                            alpha, beta, kappa, ukf_args{:});
    % =====================================================================
    
    % Show visualization:
    bennu.drawBody();
    bennu.drawLmks(visible,'MarkerSize',10);
    bennu.drawLmks_hat(X_hat(7:end,ii+1),'om','MarkerSize',10);
    grid on
    orex.draw(200,'LineWidth',2);
    view([ii*3 20])
    camva(3)
    setLims(1.1*apoapsis)
    drawnow
    
    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
%     lmks_hist(:,:,ii+1) = bennu.lmks_i;
%     visible_hist(ii,:)= visible; 
%     sig3(:,ii+1) = 3*sqrt(diag(P));
%     meas_hist(size(measurement,1),ii+1) = measurement;
end

%% Align Maps:
estimated_map = reshape(X_hat(7:end,end),[],3)';

[regParams,~, errorStats_ukf] = absor(bennu.lmks_i,estimated_map,'doScale',true,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);

figure()
scatter3(bennu.lmks_i(1,:),bennu.lmks_i(2,:),bennu.lmks_i(3,:),...
         20,'g','filled'); hold on; axis equal; grid on
scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
         100,'m','o','LineWidth',2); hold on; axis equal; grid on
disp(errorStats_ukf.errlsq)

%% Run the Smoother:


%% Plot Results:
