% clear; matlabrc; clc; close all; rng(1);
import = @(x) addpath(genpath(x));
import('data')
import('lib')
import('src')
import('filters')

%% Simulation Setup:

% Load the scenario from before:
load('results/environment_final.mat')
load('results/estimated_states.mat')
camera = orex.camera;

% Number of degrees of gravity field to estimate:
num_degrees = 4; % max degrees to estimate

% Create the ultrawideband "contellation":
num_uwb = 3;
r_uwb = repmat(orex.r,1,num_uwb);
% v_uwb = orex.v + .015*randn(3,num_uwb);
v_uwb = init_uwb_vel(r_uwb,orex.v);
    
% Extract the map and initial state estimate:
estimated_map = reshape(X_hat(7:end,end),[],3)';
num_lmks = size(estimated_map,2);
r_hat = X_hat(1:3,end);
v_hat = X_hat(4:6,end);

% Define measurement error values:
optical_std = camera.fov(1)/camera.resolution(1); %(pixels) Error in image measurements
range_uwb_std = 0.02; %(m) uncertainty in range measurements
angle_uwb_std = deg2rad(7); %(rad) uncertainty in PDoA angle measurements

% Create the ultra-wideband constellation:
uwb = UWB(r_uwb,v_uwb, range_uwb_std,angle_uwb_std);
 
% Time setup:
dt = 10*60;
orbital_period = 89836;
% duration = 3*86400;
duration = 1*orbital_period;
tspan = dt:dt:duration;


%% Filter Initialization:
alpha = 1e-3;
beta  = 2000;
kappa = 4;

% Estimation covariance initial values:
p_r_unc   = 100;
p_v_unc   = 0.1;
p_r_uwb_unc = 100;
p_v_uwb_unc = 0.1;
p_mu_unc  = 1;
p_cnm_unc = 1;
p_snm_unc = 1;

% Process noise covariance initial values:
q_r     = 1e-6;
q_v     = 1e-6;
q_r_uwb = 1e-6;
q_v_uwb = 1e-6;
q_mu    = 1e-6;
q_cnm   = 1e-6;
q_snm   = 1e-6;

% Initial estimates for all of the coefficients:
[Cnm_vec, Snm_vec] = coeffs2vec(zeros(num_degrees+1), zeros(num_degrees+1));

%% Initialize:
L = length(tspan);

r_hist = zeros(3,L); r_hist(:,1) = orex.r; 
v_hist = zeros(3,L); v_hist(:,1) = orex.v;
r_uwb_hist = zeros(3*uwb.num_transceivers,L); r_uwb_hist(:,1) = r_uwb(:);
Cnm_hist = zeros(size(Cnm_vec,1),L);
Snm_hist = zeros(size(Snm_vec,1),L);


%% Initialize Map Estimate:
% Preallocate estimate:
num_coeffs = size(Cnm_vec,1)+size(Snm_vec,1);
X_hat = zeros(6 + 6*uwb.num_transceivers + 1 + num_coeffs, L);
sig3 = X_hat;

% Store as initial estimate:
X_hat(1:6,1) = [r_hat; v_hat];
X_hat(7:6+6*uwb.num_transceivers) = [uwb.r(:) + 0*p_r_uwb_unc*randn(numel(uwb.r),1);
                                     uwb.v(:) + 0*p_v_uwb_unc*randn(numel(uwb.r),1)];
X_hat(7+6*uwb.num_transceivers,1) = bennu.mu + 0*p_mu_unc*randn;

% TODO: REMOVE THESE DEBUGGING LINES:
[Cnm_true, Snm_true] = coeffs2vec(bennu.Cnm(1:5,1:5), bennu.Snm(1:5,1:5));
X_true = [orex.r; orex.v; uwb.r(:); uwb.v(:); bennu.mu; Cnm_true; Snm_true];

% Initial Covariance matrices:
P = diag([p_r_unc*ones(3,1);
          p_v_unc*ones(3,1);
          p_r_uwb_unc*ones(3*uwb.num_transceivers,1);
          p_v_uwb_unc*ones(3*uwb.num_transceivers,1);
          p_mu_unc;
          p_cnm_unc*ones(size(Cnm_vec,1),1);
          p_snm_unc*ones(size(Snm_vec,1),1)]);
      
Q = diag([q_r*ones(3,1);
          q_v*ones(3,1);
          q_r_uwb*ones(3*uwb.num_transceivers,1);
          q_v_uwb*ones(3*uwb.num_transceivers,1);
          q_mu;
          q_cnm*ones(size(Cnm_vec,1),1);
          q_snm*ones(size(Snm_vec,1),1)]);
      
num_range_meas = 2*uwb.num_transceivers + factorial(uwb.num_transceivers)/factorial(uwb.num_transceivers - 2);
R = diag([(optical_std^2)*ones(1,2*num_lmks),...
          (optical_std^2)*ones(1,2*uwb.num_transceivers),...
          (range_uwb_std^2)*ones(1,num_range_meas)]);

sig3(:,1) = 3*sqrt(diag(P));

% Format arguments for ukf input:
dynamics_args    = {bennu, uwb.num_transceivers}; % Passes by reference
measurement_args = {orex, estimated_map, uwb.num_transceivers};
ukf_args = {dynamics_args, measurement_args};

%% Run Simulation:
% figure('units','normalized','outerposition',[0 0 1 1])
% vid = VideoWriter('results/uwb_ranges.mp4','MPEG-4');
% vid.FrameRate = 30;
% vid.Quality = 60;
% open(vid)
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    orex.rotmat_hat = orex.rotmat;
    uwb.propagate(dt,bennu);
    
    % Update the map estimate:
    estimated_map = bennu.rotate(dt)'*estimated_map;
    
    % COLLECT MEASUREMENTS: ===============================================
    % Collect measurement of the landmarks:
    [image_lmks,visible_lmks,lmk_inds] = orex.imageBody(bennu,optical_std);
    
    % Collect optical measurements of the UWB transceivers:
    [image_uwb,visible_uwb] = orex.image(uwb.r,bennu,optical_std);
    
    % Collect the UWB inter-transceiver range and angle measurements:
    [range_meas,range_avail,vec_pairs] = uwb.measureRanges(orex,bennu,range_uwb_std);
%     [angle_meas,angle_avail] = uwb.measureAngles(range_uwb_std,vec_pairs,range_avail);
    
    
    % FORMAT MEASUREMENTS FOR UKF: ========================================
    % Sort the identified landmarks by their landmark id:
    labeled_image_lmks = sortrows([bennu.lmks_lbl(visible_lmks); image_lmks]')';
    
    measurement = [labeled_image_lmks(2,:)'; labeled_image_lmks(3,:)';
                   image_uwb(1,:)'; image_uwb(2,:)';
                   range_meas];
    
    % Determine which measurements of the states are available:
    [~,ia] = intersect(1:num_lmks,bennu.lmks_lbl(visible_lmks));
    lmkAvailBools = false(num_lmks,1);
    lmkAvailBools(ia) = true;
    measAvails = [lmkAvailBools; lmkAvailBools;
                  visible_uwb(1,:)'; visible_uwb(2,:)';
                  range_avail];
    
    
    % Run the filter:======================================================
    [X_hat(:,ii+1),P] = ukf(@ukfOrbit_grav, @ukfCamera_grav,...
                            X_hat(:,ii), dt,...
                            P, Q, R, measAvails, measurement,...
                            alpha, beta, kappa, ukf_args{:});
    % =====================================================================

    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
    r_uwb_hist(:,ii+1) = uwb.r(:);
%     sig3(:,ii+1) = 3*sqrt(diag(P));
%     P_hist(1:size(P,1),1:size(P,2),ii+1) = P;
%     meas_hist(1:size(measurement,1),ii) = measurement;
%     meas_avail_hist(1:size(measAvails),ii) = measAvails;
%     pred_hist(1:size(measurement,1),ii) = ukfCamera(dt,X_hat(:,ii),orex,measAvails);

    % Only turn UWB on every 5 time steps:
    if mod(ii,1) == 0
        uwb_on = true;
    else
        uwb_on = false;
    end
    
    disp(ii/L)

    % Drawing stuff:
%     bennu.drawBody();
%     orex.draw(200,'LineWidth',2);
%     orex.drawTraj(r_hist(:,1:ii+1),'b','LineWidth',0.5);
%     uwb.draw('.m','MarkerSize',10);
%     uwb.drawTraj(r_uwb_hist(:,1:ii+1),'r','LineWidth',0.5);
%     uwb.drawConnections(orex,bennu,uwb_on,'g','LineWidth',0.5);
%     if ii == 1
%         xlim(1.3*[-2.2974, 2.7306]*1e3)
%         ylim(1.3*[-2.0187, 1.9469]*1e3)
%         zlim(1.3*[-1.1488, 2.2140]*1e3)
%     end
%     camva(4)
%     grid on
%     view([200-ii/2 20])
%     drawnow
%     frame = getframe(gcf);
%     writeVideo(vid,frame);
end
% close(vid)

%% Show the Results:
plot3(X_hat(1,:),X_hat(2,:),X_hat(3,:),'--r'); hold on
plot3(r_hist(1,:),r_hist(2,:),r_hist(3,:),'b');
axis equal
grid on

% mu_est =
% [Cnm_est, Snm_est] = vec2coeffs();

% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
%     bennu.gravityField.draw();
%     title('True 10x10 Field')
% subplot(1,2,2)
%     bennu.gravityField.draw(Cnm_est, Snm_est, mu_est);
%     title('Estimated 4x4 Field')