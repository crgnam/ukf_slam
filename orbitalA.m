clear; matlabrc; clc; close all;
import = @(x) addpath(genpath(x));
import('data')
import('lib')
import('src')

%% Simulation Setup:
% Objects to be tracked:
num_lmks = 50; % Number of features to be tracked

% Orbit setup: (currently using rough orbital A)
semi_major   = 2000; %(meters) Orbital radius
eccentricity = 0.1;
inclination  = 80;
arg_peri     = 0;
right_asc    = 0;
mean_anom    = 0;
sun_vec  = [1;0;0]; % Sun vector (for illumination visibility of lmks)
apoapsis = 2*semi_major*(1+eccentricity); 

% Setup camera (pinhole model):
f = .055; %(m) Camera focal length
resolution = [1920 1080]; %(pixels) Image dimensions
sensor     = [3 2]/1000; %(m) Sensor dimensions
fov = [sensor(1) sensor(2)]/f;
K = [f 0 0 0; 0 f 0 0; 0 0 1 0]; % Camera projection matrix:
std_meas = fov(1)/resolution(1); % Measurement uncertainty of a single pixel
 
% Time setup:
dt = 60;
duration = 86400;


%% Initialize:
% Instantiate camera and asteroid objects:
camera = Camera(K,fov,sensor,resolution);
bennu  = Asteroid('data/bennu.obj','data/bennugrav.mat',sun_vec,num_lmks,1000);

% Calculate important initial values:
tspan = dt:dt:duration;
[r,v] = orb2rv(semi_major,eccentricity,inclination,...
               arg_peri,right_asc,mean_anom,bennu.mu);
                           
% Instantiate the orex vehicle:
rotmat = nadir(r,v);
orex   = Spacecraft(r,v,rotmat,camera);

% Initialize filter:
[P,Q,R,X_hat,sig3] = ukf_init();

% Pre-allocate memory for logging data:

%% Run Simulation: