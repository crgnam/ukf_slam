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

%% Show the initial projection:
% figure('units','normalized','outerposition',[0 0 1 1])
% h1 = plot3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),'.b','MarkerSize',20); hold on
% h2 = plot3(projectedPoints(1,:),projectedPoints(2,:),projectedPoints(3,:),'xr','MarkerSize',10,'LineWidth',2);
% cam = Attitude(150,'LineWidth',2);
% cam.draw(r,rotMat)
% camva(1.5)
% legend([h1,h2],'True Map Points','Map Initial Projection','location','southeast')
% set(findall(gcf,'-property','FontSize'),'FontSize',20)
% axis equal
% grid on
% xlim([-apoapsis apoapsis])
% ylim([-apoapsis apoapsis])
% zlim([-apoapsis apoapsis])
% vid = VideoWriter('projection.mp4','MPEG-4');
% vid.Quality = 90;
% open(vid)
% for ii = 1:360
%     view([ii 20])
%     frame = getframe(gcf);
%     writeVideo(vid,frame);
% end
% close(vid)

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
    vid = VideoWriter('results/animation.mp4','MPEG-4');
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

%% Show Forward Results:  
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Kabsch algorithm to transform to truth space:
[regParams,~, errorStats_ukf] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);
M_ukf = (1/S_ukf)*R_ukf'*(X_hat(1:3,:) - t_ukf);

if save_video
    vid = VideoWriter('results/ukf_results.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'g','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             100,'m','o','LineWidth',2); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(M_ukf(1,:),M_ukf(2,:),M_ukf(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(M_ukf(1,1),M_ukf(2,1),M_ukf(3,1),'.r','MarkerSize',20)
    p = bennu.drawBodyStandalone(gca); p.FaceAlpha=0.8;
    legend('True Map','Estimated Map','True Trajectory','Estimated Trajectory',...
           'True Start','Estimated Start','Bennu','location','southeast')
    camva(2)
    xlim([-apoapsis apoapsis])
    ylim([-apoapsis apoapsis])
    zlim([-apoapsis apoapsis])
    title('Forward UKF Results')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
if save_video
    for ii = 1:360
        view([ii 20])
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    close(vid)
end

%% Plot Comparison between maps:
figure('units','normalized','outerposition',[0 0 1 1])
MS = 20;

estimated_map = reshape(X_hat(7:end,end),[],3)';
lmk_err = bennu.lmks - estimated_map;
N = size(lmk_err,2);
m = ones(1,N)/N;
lmk_lrms = sqrt(sum(sum(bsxfun(@times,m,lmk_err).*lmk_err)));
subplot(3,3,1)
    plot(lmk_err(1,:),'.b','MarkerSize',MS)
    grid on
    ylabel(['Raw (',num2str(lmk_lrms),')'])
    title('X Error (m)')
    ylim([-100 100])
    xlim([0 inf])
subplot(3,3,2)
    plot(lmk_err(2,:),'.b','MarkerSize',MS)
    grid on
    title('Y Error (m)')
    ylim([-100 100])
    xlim([0 inf])
    
subplot(3,3,3)
    plot(lmk_err(3,:),'.b','MarkerSize',MS)
    grid on
    title('Z Error')
    ylim([-100 100])
    xlim([0 inf])
    
    
estimated_map = reshape(X_hat(7:end,end),[],3)';
[R_ukf, t_ukf] = kabsch(bennu.lmks, estimated_map);
estimated_map = R_ukf'*(estimated_map - t_ukf);
lmk_err = bennu.lmks - estimated_map;
N = size(lmk_err,2);
m = ones(1,N)/N;
lmk_lrms = sqrt(sum(sum(bsxfun(@times,m,lmk_err).*lmk_err)));
subplot(3,3,4)
    plot(lmk_err(1,:),'.b','MarkerSize',MS)
    grid on
    ylabel(['Kabsch (',num2str(lmk_lrms),')'])
%     xlabel('Landmark')
    ylim([-100 100])
    xlim([0 inf])

subplot(3,3,5)
    plot(lmk_err(2,:),'.b','MarkerSize',MS)
    grid on
%     xlabel('Landmark')
    ylim([-100 100])
    xlim([0 inf])

subplot(3,3,6)
    plot(lmk_err(3,:),'.b','MarkerSize',MS)
    grid on
%     xlabel('Landmark')
    ylim([-100 100])
    xlim([0 inf])


estimated_map = reshape(X_hat(7:end,end),[],3)';
[regParams,~, errorStats_ukf] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R_ukf = regParams.R;
t_ukf = regParams.t;
S_ukf = regParams.s;
estimated_map = (1/S_ukf)*R_ukf'*(estimated_map - t_ukf);
lmk_err = bennu.lmks - estimated_map;
N = size(lmk_err,2);
m = ones(1,N)/N;
lmk_lrms = sqrt(sum(sum(bsxfun(@times,m,lmk_err).*lmk_err)));
subplot(3,3,7)
    plot(lmk_err(1,:),'.b','MarkerSize',MS)
    grid on
    ylabel(['Horn (',num2str(lmk_lrms),')'])
    xlabel('Landmark')
    ylim([-1 1])
    xlim([0 inf])

subplot(3,3,8)
    plot(lmk_err(2,:),'.b','MarkerSize',MS)
    grid on
    xlabel('Landmark')
    ylim([-1 1])
    xlim([0 inf])

subplot(3,3,9)
    plot(lmk_err(3,:),'.b','MarkerSize',MS)
    grid on
    xlabel('Landmark')
    ylim([-1 1])
    xlim([0 inf])
    
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,'results/map_compares.png')

%% Run the smoother:
[M,P_hist,D] = urts(X_hat,P_hist,@ukfOrbit,dt,Q,dynamics_args,alpha,beta,kappa,0,1);

% Get the 3-sigma bounds:
for ii = 1:size(P_hist,3)
    sig3_rts(:,ii) = 3*sqrt(diag(P_hist(:,:,ii)));
end
disp('Completed Smoother')
    
%% Show Smoother Results:
estimated_map = reshape(X_hat(7:end,end),[],3)';

% Use Horn's Method to transform to truth space:
[regParams,Bfit,errorStats_rts] = absor(bennu.lmks,estimated_map,'doScale',true,'doTrans',true);
R_rts = regParams.R;
t_rts = regParams.t;
S_rts = regParams.s;

% Transform map and states:
estimated_map = (1/S_rts)*R_rts'*(estimated_map - t_rts);
M_rts(1:3,:) = (1/S_rts)*R_rts'*(M(1:3,:) - t_rts);

if save_video
    vid = VideoWriter('results/urts_results.mp4','MPEG-4');
    vid.Quality = 90;
    open(vid)
end
figure('units','normalized','outerposition',[0 0 1 1])
    scatter3(bennu.lmks(1,:),bennu.lmks(2,:),bennu.lmks(3,:),...
             20,'g','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             100,'m','o','LineWidth',2); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(M_rts(1,:),M_rts(2,:),M_rts(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(M_rts(1,1),M_rts(2,1),M_rts(3,1),'.r','MarkerSize',20)
    p = bennu.drawBodyStandalone(gca); p.FaceAlpha=0.8;
    legend('True Map','Estimated Map','True Trajectory','Estimated Trajectory',...
           'True Start','Estimated Start','Bennu','location','southeast')
    camva(2)
    xlim([-apoapsis apoapsis])
    ylim([-apoapsis apoapsis])
    zlim([-apoapsis apoapsis])
    title('Unscented RTS Smoother Results')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    if save_video
        for ii = 1:360
            view([ii 20])
            frame = getframe(gcf);
            writeVideo(vid,frame);
        end
        close(vid)
    end
        
%% Calculate Residuals:
pt_resid = nan(size(meas_history));
pf_resid = nan(size(meas_history));
std_pt_row = zeros(1,size(pf_resid,2));
std_pt_line = std_pt_row;
std_pf_row = std_pt_row;
std_pf_line = std_pt_row;
for ii = 1:size(X_hat,2)
    % Get the measurement availability:
    avails = avails_history(:,ii)==1;
    
    % Calculate pass-through predicted measurement:
    pt_predict = ukfCamera(dt, X_hat(:,ii),rotMat_hist(:,:,ii),K, avails);
    pf_predict = ukfCamera(dt, M(:,ii),    rotMat_hist(:,:,ii),K, avails);
    
    % Calculate the measurement residuals:
    nn = sum(avails);
    if nn ~= 0
        pt_resid(1:nn/2,ii) = pt_predict(1:nn/2) - meas_history(1:nn/2,ii);
        pf_resid(1:nn/2,ii) = pf_predict(1:nn/2) - meas_history(1:nn/2,ii);
        pt_resid(N+1:N+nn/2,ii) = pt_predict(nn/2+1:end) - meas_history(nn/2+1:nn,ii);
        pf_resid(N+1:N+nn/2,ii) = pf_predict(nn/2+1:end) - meas_history(nn/2+1:nn,ii);

        % Calculate statistics:
        std_pt_row(ii)  = std(pt_resid(1:nn/2,ii));
        std_pt_line(ii) = std(pt_resid(N+1:N+nn/2,ii));
        std_pf_row(ii)  = std(pf_resid(1:nn/2,ii));
        std_pf_line(ii) = std(pf_resid(N+1:N+nn/2,ii));
    end
end

% Plot the residuals:
pt_line = pt_resid(1:N,:);
pt_row  = pt_resid(N+1:end,:);
pf_line = pf_resid(1:N,:);
pf_row  = pf_resid(N+1:end,:);

%% Plot Residuals:
MS = 1;
YLIM  = 2;
YLIM2 = 1;
NN = 5;
figure('units','normalized','outerposition',[0 0 1 1])
    t_plt = repmat((0:dt:tspan(end))/60,N,1);
    subplot(2,2,1)
        plot(t_plt,pt_line*res_x,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_line*res_x,NN)
        ylim([-YLIM YLIM])
        xlim([0 inf])
        title('Forward UKF')
        ylabel('Line Error (pixels)')
        grid on
    subplot(2,2,2)
        plot(t_plt,pf_line*res_x,'.m','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pf_line*res_x,NN)
        ylim([-YLIM YLIM])
        xlim([0 inf])
        title('RTS Smoother')
        grid on

    subplot(2,2,3)
        plot(t_plt,pt_row*res_y,'.b','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pt_row*res_y,NN)
        ylim([-YLIM2 YLIM2])
        xlim([0 inf])
        ylabel('Row Error (pixels)')
        xlabel('Time (min)')
        grid on
    subplot(2,2,4)
        plot(t_plt,pf_row*res_y,'.m','MarkerSize',MS); hold on
        drawBounds(t_plt,std_pf_row*res_y,NN)
        xlabel('Time (min)')
        ylim([-YLIM2 YLIM2])
        xlim([0 inf])
        grid on
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,'results/residuals.png')

%% Calculate Trajectory Residuals:
% Get missing parts:
Mv_ukf = (1/S_ukf)*R_ukf'*(X_hat(4:6,:) - t_ukf);
r_ukf_resid = M_ukf - X(1:3,:);
v_ukf_resid = Mv_ukf - X(4:6,:);
r_sig_ukf = (1/3)*(1/S_ukf)*sig3(1:3,:);
v_sig_ukf = (1/3)*(1/S_ukf)*sig3(4:6,:);

Mv_rts = (1/S_rts)*R_rts'*(M(4:6,:) - t_rts);
r_rts_resid = M_rts  - X(1:3,:);
v_rts_resid = Mv_rts - X(4:6,:);
r_sig_rts = (1/3)*(1/S_rts)*sig3_rts(1:3,:);
v_sig_rts = (1/3)*(1/S_rts)*sig3_rts(4:6,:);

%% Plot Trajectory Residuals:
figure('units','normalized','outerposition',[0 0 1 1])
LW = 2;
NN = 5;
subplot(3,2,1)
    plot(t_plt,r_ukf_resid(1,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_ukf(1,:),NN)
    grid on
    ylabel('r_x Error (m)')
    title('Forward UKF')
    xlim([0 inf])
    
subplot(3,2,2)
    plot(t_plt,r_rts_resid(1,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(1,:),NN)
    grid on
    title('RTS Smoother')
    xlim([0 inf])
    
subplot(3,2,3)
    plot(t_plt,r_ukf_resid(2,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_ukf(2,:),NN)
    grid on
    ylabel('r_y Error (m)')
    xlim([0 inf])
    
subplot(3,2,4)
    plot(t_plt,r_rts_resid(2,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(2,:),NN)
    grid on
    xlim([0 inf])
    
subplot(3,2,5)
    plot(t_plt,r_ukf_resid(3,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_ukf(3,:),NN)
    grid on
    ylabel('r_z Error (m)')
    xlabel('Time (min)')
    xlim([0 inf])
    
subplot(3,2,6)
    plot(t_plt,r_rts_resid(3,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,r_sig_rts(3,:),NN)
    grid on
    xlabel('Time (min)')
    xlim([0 inf])
    
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,'results/traj_error.png')

%% Plot Velocity Residuals:
figure('units','normalized','outerposition',[0 0 1 1])
LW = 2;
NN = 5;
subplot(3,2,1)
    plot(t_plt,v_ukf_resid(1,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_ukf(1,:),NN)
    grid on
    ylabel('v_x Error (m/s)')
    title('Forward UKF')
    xlim([0 inf])
    
subplot(3,2,2)
    plot(t_plt,v_rts_resid(1,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_rts(1,:),NN)
    grid on
    title('RTS Smoother')
    xlim([0 inf])
    
subplot(3,2,3)
    plot(t_plt,v_ukf_resid(2,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_ukf(2,:),NN)
    grid on
    ylabel('v_y Error (m/s)')
    xlim([0 inf])
    
subplot(3,2,4)
    plot(t_plt,v_rts_resid(2,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_rts(2,:),NN)
    grid on
    xlim([0 inf])
    
subplot(3,2,5)
    plot(t_plt,v_ukf_resid(3,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_ukf(3,:),NN)
    grid on
    ylabel('v_z Error (m/s)')
    xlabel('Time (min)')
    xlim([0 inf])
    
subplot(3,2,6)
    plot(t_plt,v_rts_resid(3,:),'b','LineWidth',LW); hold on
    drawBounds(t_plt,v_sig_rts(3,:),NN)
    grid on
    xlabel('Time (min)')
    xlim([0 inf])
    
set(findall(gcf,'-property','FontSize'),'FontSize',20)
saveas(gcf,'results/vel_error.png')

%% Plot errors