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
duration = 6*86400;

%% Filter Initialization:
alpha = 1e-5;
beta  = 2000;
kappa = 4;

% Estimation covariance initial values:
p_r_unc = 100;
p_v_unc = 1;
p_lmk_unc = [100; 100; 100]; % Uncertainties in lmks x,y,z

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
orex.rotmat_hat = nadir(r_hat,v_hat);
radius_hat      = camera.radius_estimate(r_hat,orex.r,bennu.radius);
projectedPoints = initializeMap(camera,image_lmks,orex.rotmat,r_hat,radius_hat);

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
meas_hist = zeros(2*num_lmks,L);
pred_hist = zeros(2*num_lmks,L);
detection_hist = nan(1,num_lmks); dx = 1;
meas_avail_hist = zeros(2*num_lmks,L);
save('results/original.mat','orex','bennu')
for ii = 1:L-1
    % Move the scene forward one time step:
    bennu.update(dt);
    orex.propagate(dt,bennu);
    orex.rotmat_hat = nadir(X_hat(1:3,ii),X_hat(4:6,ii));
    
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
    
    % Log data:
    r_hist(:,ii+1) = orex.r;
    v_hist(:,ii+1) = orex.v;
    sig3(:,ii+1) = 3*sqrt(diag(P));
    P_hist(1:size(P,1),1:size(P,2),ii+1) = P;
    meas_hist(1:size(measurement,1),ii) = measurement;
    meas_avail_hist(1:size(measAvails),ii) = measAvails;
    pred_hist(1:size(measurement,1),ii) = ukfCamera(dt,X_hat(:,ii),orex,measAvails);
    if new_detection
        detection_hist(dx) = ii;
        dx=dx+1;
    end
end
detection_hist(isnan(detection_hist)) = [];
disp('MAP DISCOVERY COMPLETE')
save('results/environment.mat','bennu','orex')

% Extract the information for the final portion of the orbit:
tspan2 = dt:dt:86400;
L2 = length(tspan2);
X_hat2 = X_hat(:,end-L2:end);
sig32  = sig3(:,end-L2:end);
meas_hist2 = meas_hist(:,end-L2:end);
P_hist2 = P_hist(:,:,end-L2:end);
r_hist2 = r_hist(:,end-L2:end);
meas_avail_hist2 = meas_avail_hist(:,end-L2:end);

%% Map Alignment for Forward UKF:
% load('results/environment.mat')
% % Format the estimated map from the UKF output:
% estimated_map = reshape(X_hat(7:end,end),[],3)';
% 
% % Get the corresponding truth map from the identified/labelled Bennu
% % features:
% true_map_lbl = sortrows([bennu.lmks_lbl(bennu.lmks_obs); bennu.lmks_i(:,bennu.lmks_obs)]')';
% true_map = true_map_lbl(2:end,:);
% 
% % Use Horn's algorithm to perform map alignment:
% [regParams,~, errorStats_ukf] = absor(true_map,estimated_map,'doScale',true,'doTrans',true);
% R_ukf = regParams.R;
% t_ukf = regParams.t;
% S_ukf = regParams.s;
% estimated_map2 = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);
% traj = (1/S_ukf)*R_ukf'*(X_hat(1:3,:) - t_ukf);

% figure('units','normalized','outerposition',[0 0 1 1])
% C = [linspace(1,0,size(true_map,2))',linspace(0,0,size(true_map,2))' ,linspace(0,1,size(true_map,2))'];
% scatter3(true_map(1,:),true_map(2,:),true_map(3,:),...
%          20,'g','filled'); hold on; axis equal; grid on
% scatter3(estimated_map2(1,:),estimated_map2(2,:),estimated_map2(3,:),...
%          100,'m','o','LineWidth',2); hold on; axis equal; grid on
% plot3(traj(1,:),traj(2,:),traj(3,:),'r')
% plot3(traj(1,1),traj(2,1),traj(3,1),'.r','MarkerSize',20)
% plot3(traj(1,end),traj(2,end),traj(3,end),'xr','MarkerSize',20,'LineWidth',2)
% plot3(r_hist2(1,:),r_hist2(2,:),r_hist2(3,:),'b')
% camva(3)
% setLims(2.5*apoapsis)
% drawnow
% disp(errorStats_ukf.errlsq)
% disp(errorStats_ukf.errmax)

%% Show the UKF 3-sigma bounds:
% figure()
% LW = 2;
% NN = 5;
% r_ukf_resid = X_hat(1:3,:) - r_hist;
% y1 = 50;
% y2 = 20;
%     
% t_plt = (dt:dt:tspan(end))/60;
% subplot(3,1,1)
%     plot(t_plt,r_ukf_resid(1,:),'b','LineWidth',LW); hold on
%     drawBounds(t_plt,sig3(1,:),NN)
%     grid on
%     title('Forward UKF')
%     ylabel('error_x (m)')
%     ylim([-y1 y1])
%     xlim([0 inf])
%     
% subplot(3,1,2)
%     plot(t_plt,r_ukf_resid(2,:),'b','LineWidth',LW); hold on
%     drawBounds(t_plt,sig3(2,:),NN)
%     grid on
%     xlim([0 inf])
%     ylabel('error_y (m)')
%     ylim([-y2 y2])
% 
% subplot(3,1,3)
%     plot(t_plt,r_ukf_resid(3,:),'b','LineWidth',LW); hold on
%     drawBounds(t_plt,sig3(3,:),NN)
%     grid on
%     xlabel('Time (min)')
%     ylabel('error_z (m)')
%     xlim([0 inf])
%     ylim([-y1 y1])
%     
% set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Show the Measurement Residuals:
pt_resid = nan(size(meas_hist));
std_pt_row = zeros(1,size(pt_resid,2));
std_pt_line = std_pt_row;
for ii = 1:size(X_hat,2)
    % Get the measurement availability:
    nn = nnz(meas_hist(:,ii));
    
    % Calculate the measurement residuals:
    if nn ~= 0
        pt_resid(1:nn/2,ii)     = pred_hist(1:nn/2,ii) - meas_hist(1:nn/2,ii);
        pt_resid(num_lmks+1:num_lmks+nn/2,ii) = pred_hist(nn/2+1:nn,ii) - meas_hist(nn/2+1:nn,ii);

        % Calculate statistics:
        std_pt_row(ii)  = std(pt_resid(1:nn/2,ii));
        std_pt_line(ii) = std(pt_resid(num_lmks+1:num_lmks+nn/2,ii));
    end
end

% Extract and format information for plotting:
pt_line = pt_resid(1:num_lmks,:);
pt_row  = pt_resid(num_lmks+1:end,:);


%% Plot the residuals:
MS = 5;
YLIM  = 2;
YLIM2 = 2;
NN = 5;

res_x = resolution(1);
res_y = resolution(2);
figure('units','normalized','outerposition',[0 0 1 1])
    t_plt = repmat((dt:dt:tspan(end))/60,num_lmks,1);
    subplot(2,1,1)
        plot(t_plt,pt_line*res_x,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_line*res_x,NN)
        ylim([-YLIM YLIM])
        drawDetections(t_plt,detection_hist,ylim,'--r','LineWidth',1)
        xlim([0 inf])
        title('Forward UKF')
        ylabel('Line Error (pixels)')
        grid on
        grid on

    subplot(2,1,2)
        plot(t_plt,pt_row*res_y,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_row*res_y,NN)
        ylim([-YLIM2 YLIM2])
        drawDetections(t_plt,detection_hist,ylim,'--r','LineWidth',1)
        xlim([0 inf])
        ylabel('Row Error (pixels)')
        xlabel('Time (min)')
        grid on
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,'results/residuals.png')

%% Run the Smoother:
disp('Running Smoother')
load('results/environment.mat')
[M,P_hist,D] = urts(X_hat2,P_hist2,@ukfOrbit_RTS,dt,Q,dynamics_args,alpha,beta,kappa,0,1);

r_sig_rts = zeros(size(P_hist,1),size(P_hist,3));
for ii = 1:size(P_hist,3)
   r_sig_rts(:,ii) = 3*sqrt(diag(P_hist(:,:,ii)));
end

%% Generate predicted measurements from Smoother states:
disp('Generating Predicted Measurements')
load('results/original.mat')
% Move scenario forward to the starting location:
L_adjust = length(tspan(1:end-L2));
for ii = 1:L_adjust-2
    bennu.update(dt);
    orex.propagate(dt,bennu);
end

pred_hist2 = zeros(size(meas_hist2));
for ii = 1:size(X_hat2,2)    
    bennu.update(dt);
    orex.propagate(dt,bennu);
    
    nn = nnz(meas_hist2(:,ii));
    measAvails = meas_avail_hist2(:,ii)==1;
    pred_hist2(1:nn,ii) = ukfCamera(dt,M(:,ii),orex,measAvails);
%     plot(pred_hist2(1:nn/2,ii),pred_hist2(1+nn/2:nn,ii),'xr'); hold on
%     plot(meas_hist2(1:nn/2,ii),meas_hist2(1+nn/2:nn,ii),'.b');
end

% Get the post-fit residuals:
pt_resid = nan(size(meas_hist2));
std_pt_row = zeros(1,size(pt_resid,2));
std_pt_line = std_pt_row;
for ii = 1:size(X_hat2,2)
    % Get the measurement availability:
    nn = nnz(meas_hist2(:,ii));
    
    % Calculate the measurement residuals:
    if nn ~= 0
        pt_resid(1:nn/2,ii)     = pred_hist2(1:nn/2,ii) - meas_hist2(1:nn/2,ii);
        pt_resid(num_lmks+1:num_lmks+nn/2,ii) = pred_hist2(nn/2+1:nn,ii) - meas_hist2(nn/2+1:nn,ii);

        % Calculate statistics:
        std_pt_row(ii)  = std(pt_resid(1:nn/2,ii));
        std_pt_line(ii) = std(pt_resid(num_lmks+1:num_lmks+nn/2,ii));
    end
end

% Extract and format information for plotting:
pt_line = pt_resid(1:num_lmks,:);
pt_row  = pt_resid(num_lmks+1:end,:);

% Show the measurement residuals:
MS = 4;
YLIM  = 2;
YLIM2 = 2;
NN = 5;

res_x = resolution(1);
res_y = resolution(2);
figure('units','normalized','outerposition',[0 0 1 1])
    t_plt = repmat((0:dt:tspan2(end))/60,num_lmks,1);
    subplot(2,1,1)
        plot(t_plt,pt_line*res_x,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_line*res_x,NN)
        ylim([-YLIM YLIM])
        xlim([0 inf])
        title('URTS SMoother')
        ylabel('Line Error (pixels)')
        grid on
        grid on

    subplot(2,1,2)
        plot(t_plt,pt_row*res_y,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_row*res_y,NN)
        ylim([-YLIM2 YLIM2])
        xlim([0 inf])
        ylabel('Row Error (pixels)')
        xlabel('Time (min)')
        grid on
set(findall(gcf,'-property','FontSize'),'FontSize',20)


%% Align Maps:
load('results/environment.mat')
% Format the estimated map from the UKF output:
estimated_map = reshape(M(7:end,end),[],3)';

% Get the corresponding truth map from the identified/labelled Bennu
% features:
true_map_lbl = sortrows([bennu.lmks_lbl(bennu.lmks_obs); bennu.lmks_i(:,bennu.lmks_obs)]')';
true_map = true_map_lbl(2:end,:);

% Use Horn's algorithm to perform map alignment:
[regParams,~, errorStats_ukf] = absor(true_map,estimated_map,'doScale',false,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map2 = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);
traj = (1/S_ukf)*R_ukf'*(M(1:3,:) - t_ukf);

% figure('units','normalized','outerposition',[0 0 1 1])
% C = [linspace(1,0,size(true_map,2))',linspace(0,0,size(true_map,2))' ,linspace(0,1,size(true_map,2))'];
% scatter3(true_map(1,:),true_map(2,:),true_map(3,:),...
%          20,'g','filled'); hold on; axis equal; grid on
% scatter3(estimated_map2(1,:),estimated_map2(2,:),estimated_map2(3,:),...
%          100,'m','o','LineWidth',2); hold on; axis equal; grid on
% plot3(traj(1,:),traj(2,:),traj(3,:),'r')
% plot3(traj(1,1),traj(2,1),traj(3,1),'.r','MarkerSize',20)
% plot3(traj(1,end),traj(2,end),traj(3,end),'xr','MarkerSize',20,'LineWidth',2)
% plot3(r_hist2(1,:),r_hist2(2,:),r_hist2(3,:),'b')
% camva(3)
% setLims(2.5*apoapsis)
% drawnow
disp(errorStats_ukf.errlsq)
disp(errorStats_ukf.errmax)


%% Show Results:
figure()
LW = 2;
NN = 5;
% r_rts_resid = M(1:3,:) - r_hist2;
r_rts_resid = traj - r_hist2;
y1 = 50;
y2 = 20;
    
t_plt = (0:dt:tspan2(end))/60;
subplot(3,1,1)
    plot(t_plt,r_rts_resid(1,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(1,:),NN)
    grid on
    title('RTS Smoother')
    ylabel('error_x (m)')
    ylim([-y1 y1])
    xlim([0 inf])
    
subplot(3,1,2)
    plot(t_plt,r_rts_resid(2,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(2,:),NN)
    grid on
    xlim([0 inf])
    ylabel('error_y (m)')
    ylim([-y2 y2])

subplot(3,1,3)
    plot(t_plt,r_rts_resid(3,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(3,:),NN)
    grid on
    xlabel('Time (min)')
    ylabel('error_z (m)')
    xlim([0 inf])
    ylim([-y1 y1])
    
set(findall(gcf,'-property','FontSize'),'FontSize',20)
