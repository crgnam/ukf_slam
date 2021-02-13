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
bennu_w             = 2*pi/(4.3*3600); %(rad/s) Rotation rate of bennu

% Setup camera (pinhole model):
f = .055; %(m) Camera focal length
resolution = [1920 1080]; %(pixels) Image dimensions
sensor     = [3 2]/1000; %(m) Sensor dimensions
fov = [sensor(1) sensor(2)]/f;
K = [f 0 0 0; 0 f 0 0; 0 0 1 0]; % Camera projection matrix:
std_meas = fov(1)/resolution(1); % Measurement uncertainty of a single pixel
 
% Time setup:
dt = 0.5*60;
duration = 1*86400;

%% Filter Initialization:
alpha = 1e-5;
beta  = 3000;
kappa = 4-156;

r_unc = 10;
v_unc = 10;
P = diag([r_unc*ones(1,3),...
          v_unc*ones(1,3),...
          50*ones(1,3*num_lmks)]);

Q = diag([1e-9*ones(1,3),...
          1e-6*ones(1,3),...
          1e-6*ones(1,3*num_lmks)]);

R = diag((std_meas^2)*ones(1,2*num_lmks));


%% Initialize:
% Instantiate camera and asteroid objects:
camera = Camera(K,fov,sensor,resolution);
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

X_hat = zeros(6 + 3*num_lmks, length(tspan));
sig3 = X_hat;
sig3_rts = sig3;
sig3(:,1) = 3*sqrt(diag(P));
P_hist = zeros(size(P,1),size(P,2),length(tspan));

%% Initialize Map Estimate:
r_hat = r + r_unc*randn(3,1);
v_hat = v + v_unc*randn(3,1);
% Take first image:
[image_lmks,visible] = orex.image(bennu,0);

% Get initial map: TEMPORARY BYPASS BECAUSE IM LAZY
% TODO: REMOVE THIS!!!!!
projectedPoints = zeros(size(bennu.lmks_i));
projectedPoints(:,visible) = bennu.lmks_i(:,visible);

% Store as initial estimate:
X_hat(:,1) = [r; v;
              projectedPoints(1,:)';
              projectedPoints(2,:)';
              projectedPoints(3,:)'];

% Format arguments for ukf input:
dynamics_args    = {bennu};
measurement_args = {orex};
ukf_args = {dynamics_args, measurement_args};

%% Run Simulation:
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    orex.rotmat = nadir(orex.r,orex.v);
    
    % Collect measurement:
    [image_lmks,visible] = orex.image(bennu,std_meas);
    
    % Run the filter:======================================================
    measurement = [image_lmks(1,:)'; image_lmks(2,:)'];
    measAvails  = [visible visible]';    
    [X_hat(:,ii+1), P] = ukf(@ukfOrbit, @ukfCamera,...
                             X_hat(:,ii), dt,...
                             P, Q, R, measAvails, measurement,...
                             alpha, beta, kappa, ukf_args{:});
    % =====================================================================
    
    % Show visualization:
%     bennu.drawBody();
    bennu.drawLmks(visible,'MarkerSize',10);
    bennu.drawLmks_hat(X_hat(7:end,ii+1),'om','MarkerSize',10)
    grid on
%     orex.draw(200,'LineWidth',2);
%     view([ii/3 20])
%     setLims(1.1*apoapsis)
    drawnow
    
    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
    lmks_hist(:,:,ii+1) = bennu.lmks_i;
    visible_hist(ii,:)= visible; 
    sig3(:,ii+1) = 3*sqrt(diag(P));
    P_hist(:,:,ii+1) = P;
    avails_hist(:,ii+1) = measAvails;
    meas_hist(1:sum(measAvails),ii+1) = measurement;
end

%% Align Maps:
estimated_map = reshape(X_hat(7:end,end),[],3)';

[regParams,~, errorStats_ukf] = absor(bennu.lmks_i,estimated_map,'doScale',true,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);

scatter3(bennu.lmks_i(1,:),bennu.lmks_i(2,:),bennu.lmks_i(3,:),...
         20,'g','filled'); hold on; axis equal; grid on
scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
         100,'m','o','LineWidth',2); hold on; axis equal; grid on
disp(errorStats_ukf.errlsq)

%% Run the Smoother:


%% Plot Results:
