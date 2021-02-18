clear; matlabrc; clc; close all; rng(1);
import = @(x) addpath(genpath(x));
import('data')
import('lib')
import('src')
import('filters')

%% Simulation Setup:

% Load the scenario from before:
load('results/environment_initial.mat')
load('results/estimated_states.mat')
 
% Time setup:
dt = 1*60;
orbital_period = 89836;
% duration = 3*86400;
duration = orbital_period;
tspan = dt:dt:duration;

%% Filter Initialization:
alpha = 1e-5;
beta  = 2000;
kappa = 4;

% Estimation covariance initial values:
p_r_unc = 10;
p_v_unc = 0.1;

% Process noise covariance initial values:
q_r = 1e-12;
q_v = 1e-12;


%% Pre-allocate memory for logging data:
L = length(tspan);

r_hist = zeros(3,L); r_hist(:,1) = r; 
v_hist = zeros(3,L); v_hist(:,1) = v;

%% Initialize Map Estimate:
r_hat = orex.r + p_r_unc*randn(3,1);
v_hat = orex.v + p_v_unc*randn(3,1);

% Preallocate estimate:
X_hat = zeros(6+3*num_tracking,L);
P_hist = zeros(6+3*num_lmks,6+3*num_lmks,L);
sig3 = X_hat;

% Store as initial estimate:
X_hat(:,1) = [r_hat; v_hat; projectedPoints];

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

%% Run Simulation:
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    orex.rotmat_hat = orex.rotmat;
    
    % Collect measurement:
    [image_lmks,visible,lmk_inds] = orex.image(bennu,std_meas);
    
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

    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
    sig3(:,ii+1) = 3*sqrt(diag(P));
    P_hist(1:size(P,1),1:size(P,2),ii+1) = P;
    meas_hist(1:size(measurement,1),ii) = measurement;
    meas_avail_hist(1:size(measAvails),ii) = measAvails;
    pred_hist(1:size(measurement,1),ii) = ukfCamera(dt,X_hat(:,ii),orex,measAvails);
end