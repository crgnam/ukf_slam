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
bennu_w             = 1*2*pi/(4.3*3600); %(rad/s) Rotation rate of bennu

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
duration = 5*86400;

%% Filter Initialization:
alpha = 1e-5;
beta  = 2000;
kappa = 4;

% Estimation covariance initial values:
p_r_unc = 50;
p_v_unc = 1;
p_lmk_unc = [250; 250; 250]; % Uncertainties in lmks x,y,z

% Process noise covariance initial values:
q_r = 1e-12;
q_v = 1e-12;
q_lmk = [1e-12; 1e-12; 1e-12];


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

X_hat(:,1) = [r_hat; v_hat; projectedPoints];

figure('units','normalized','outerposition',[0 0 1 1])
orex.drawRays(image_lmks,radius_hat,r_hat);
bennu.drawBody();
set(bennu.ptch,'AmbientStrength',.8,'FaceAlpha',0.5)
bennu.drawLmks(visible,'MarkerSize',30);
bennu.drawLmks_hat(X_hat(7:end,1),'ob','MarkerSize',10);
setLims(apoapsis)
camva(2)