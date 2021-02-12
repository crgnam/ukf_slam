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
eccentricity = 0.2;
inclination  = 80;
arg_peri     = 0;
right_asc    = 0;
mean_anom    = 0;
apoapsis = semi_major*(1+eccentricity); 

% Bennu states:
sun_vec           = [0;1;0]; % Sun vector (for illumination visibility of lmks)
bennu_orientation = ea2rotmat(5,3,0,'321'); % Orientation of Bennu
bennu_w           = 2*pi/(4.3*3600); %(rad/s) Rotation rate of bennu

% Setup camera (pinhole model):
f = .055; %(m) Camera focal length
resolution = [1920 1080]; %(pixels) Image dimensions
sensor     = [3 2]/1000; %(m) Sensor dimensions
fov = [sensor(1) sensor(2)]/f;
K = [f 0 0 0; 0 f 0 0; 0 0 1 0]; % Camera projection matrix:
std_meas = fov(1)/resolution(1); % Measurement uncertainty of a single pixel
 
% Time setup:
dt = 60;
duration = 3*86400;

% Filter Initialization:



%% Initialize:
% Instantiate camera and asteroid objects:
camera = Camera(K,fov,sensor,resolution);
bennu  = Asteroid('data/bennu.obj','data/bennugrav.mat',...
                  bennu_orientation,bennu_w,sun_vec,num_lmks,1000);

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


%% Run Simulation:
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    
    % Set nadir pointing attitude:
    orex.rotmat = nadir(orex.r,orex.v);
    
    % Collect measurement:
%     [lmks,visible] = orex.image(bennu);
    
    % Run the filter:======================================================
    
    
    % =====================================================================
    
    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
end

%% Run the Smoother:


%% Plot Results:
plot3(r_hist(1,:),r_hist(2,:),r_hist(3,:)); hold on; axis equal