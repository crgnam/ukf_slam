%% Case 2: Estimating Rotating Body, No Measurement Noise
% In this case, we also aim to estimate the rotation vector and rotation
% rate of the body.  It does this through estimating a set of euler angles
% that represent the rotation state of the asteroid
clear; matlabrc; clc; close all
addpath(genpath('src'))
addpath(genpath('ukf'))
rng(1)

% Setup:
N = 30; % Number of features to track
dt = 60;
tspan = dt:dt:4*86400;
w_bennu = 2*pi/(4.3*60*60); %(rad/s) Angular rate of Bennu
bennu_theta = 2.32;
bennu_phi = 3.12;
bennu_psi = 0;

f = .055; %(m) Camera focal length
sx = 3; %sensor x in mm
sy = 2; %sensor y in mm
fov = [sx/1000 sy/1000]/f;
K = [f 0 0 0;
     0 f 0 0;
     0 0 1 0];

% Generate fake map:
az = 2*pi*rand(1,N);
el = -pi/2 + pi*rand(1,N);
r = 250*ones(1,N) + 20*randn(1,N); %(meters) Radius of Bennu
[x,y,z] = sph2cart(az,el,r);
colors = rand(length(x),3); % Colors of each feature (unique)
% Colors (blue->red based on latitude)
colors(:,1) = 1 - (z+25)/500;
colors(:,2) = 0;
colors(:,3) = (z+250)/500;

% Setup orbital parameters:
G = 6.67430*10^-11; %(m3⋅kg–1⋅s–2) Gravitational Constant
M = (7.329*10^10); %(kg) Mass of Bennu
mu = G*M; % Standard gravitational constant of Bennu
a = 700; %(meters) Orbital radius
e = 0.3;
inc = 45;
w = 0; RAAN = 0; M0 = 0;
[r, v] = orbitalElements2PosVel(a,e,inc,w,RAAN,M0,mu);
apo = a*(1+e); 

% Preallocate memory:
X = zeros(6,length(tspan));
X(:,1) = [r; v];

% Initialize UKF:
P = diag([10*ones(1,3),...
          0.1*ones(1,3),...
          1,1,1,.05,...
          250*ones(1,3*N)]); % Initial State Covariance Matrix
Q = diag([1e-9*ones(1,3),...
          1e-9*ones(1,3),...
          1e-9*ones(1,4),...
          1e-9*ones(1,3*N)]); % Process Noise Covariance Matrix

% TODO: Hardocded measurement covariance matrix, assuming ALL features are
% visible at every time step (valid for first pass)
pix_sig = 0.001;
R = diag((pix_sig^2)*ones(1,2*N));

% Initial estimate of the map:
rotMat = nadir(r,v);
imagePoints = camera(rotMat,r,K,fov, [x;y;z]);
rays = generateRays(imagePoints,rotMat,K);
projectedPoints = r + norm(r)*rays;

% Initialize estimates:
theta_hat = 2;
phi_hat   = 3;
psi_hat   = 0;
w_hat     = 2*pi/(4*60*60);
bennu_psi = zeros(1,length(tspan));
X_hat = zeros(6 + 4 + 3*N, length(tspan));
X_hat(:,1) = [r; v;
              theta_hat; phi_hat; psi_hat; w_hat;
              projectedPoints(1,:)';
              projectedPoints(2,:)';
              projectedPoints(3,:)'];
sig3 = zeros(size(X_hat));
sig3(:,1) = 3*sqrt(diag(P));
          
% Initialize UKF model arguments:
dynamics_args = {mu};

% UKF tuning parameters:
alpha = 1;
beta  = 10;
kappa = 2;

%% Run Filter:
% For the purposes here of this test, it is assumed that ALL features are
% ALWAYS visible
for ii = 1:length(tspan)
    % Simulate the orbit:
    X(:,ii+1) = rk4(@orbitalDynamics,dt,X(:,ii),mu);
    r = X(1:3,ii+1);
    v = X(4:6,ii+1);
    
    % Define a nadir pointing attitude:
    rotMat = nadir(r,v);
    
    % Rotate the landmarks on Bennu's surface:
    bennu_psi(ii+1) = w_bennu*ii*dt;
    bennu_rotMat = ea2rotmat123(bennu_theta,bennu_phi,bennu_psi(ii+1));
    lmks = bennu_rotMat*[x;y;z];
    
    % Calculate the measurement:
    imagePoints = camera(rotMat,r,K,fov, lmks); %(CURRENTLY NO NOISE)
    measurement = [imagePoints(1,:)'; imagePoints(2,:)'];
    
    % Formulate UKF inputs:
    measAvails = true(size(measurement));
    measurement_args = {rotMat,K};
    ukf_args = {dynamics_args, measurement_args};
    
    % Run the UKF:
    [X_hat(:,ii+1), P] = ukf(@ukfOrbit_bennu_rotation, @ukfCamera_bennu_rotation,...
                             X_hat(:,ii), dt,...
                             P, Q, R, measAvails, measurement,...
                             alpha, beta, kappa, ukf_args{:});
    sig3(:,ii+1) = 3*sqrt(diag(P));
end

%% Show the Final Estimated Map:   
estimated_map = reshape(X_hat(11:end,end),[],3)';

figure()
    scatter3(x,y,z,30,'b','filled'); hold on; axis equal; grid on
    scatter3(estimated_map(1,:),estimated_map(2,:),estimated_map(3,:),...
             30,'r','x'); hold on; axis equal; grid on
    plot3(X(1,:),X(2,:),X(3,:),'b')
    plot3(X_hat(1,:),X_hat(2,:),X_hat(3,:),'r')
    plot3(X(1,1),X(2,1),X(3,1),'.b','MarkerSize',20)
    plot3(X_hat(1,1),X_hat(2,1),X_hat(3,1),'.r','MarkerSize',20)
    legend('truth map','estimated map','truth trajectory','estimated trajectory')
    
    scale = 300;
    R = scale*ea2rotmat123(bennu_theta, bennu_phi, bennu_psi(end));
    plot3([0 R(1,1)],[0 R(1,2)],[0 R(1,3)],'r');
    plot3([0 R(2,1)],[0 R(2,2)],[0 R(2,3)],'g');
    plot3([0 R(3,1)],[0 R(3,2)],[0 R(3,3)],'b');
    
    R = scale*ea2rotmat123(X_hat(7,end), X_hat(8,end), X_hat(9,end));
    plot3([0 R(1,1)],[0 R(1,2)],[0 R(1,3)],'--r');
    plot3([0 R(2,1)],[0 R(2,2)],[0 R(2,3)],'--g');
    plot3([0 R(3,1)],[0 R(3,2)],[0 R(3,3)],'--b');
    
figure()
subplot(3,1,1)
    plot(X_hat(7,:)-bennu_theta); hold on
    plot(-sig3(7,:),'--k')
    plot(sig3(7,:),'--k')
    ylabel('theta')
subplot(3,1,2)
    plot(X_hat(8,:)-bennu_phi); hold on
    plot(-sig3(8,:),'--k')
    plot(sig3(8,:),'--k')
    ylabel('phi')
% subplot(2,2,3)
%     plot(X_hat(9,:)-bennu_psi); hold on
%     plot(-sig3(9,:),'--k')
%     plot(sig3(9,:),'--k')
subplot(3,1,3)
    plot(X_hat(10,:)-w_bennu); hold on
    plot(-sig3(10,:),'--k')
    plot(sig3(10,:),'--k')
    ylim([-2*mean(sig3(10,:)) 2*mean(sig3(10,:))])
    ylabel('\omega')

% Plot the post fit residuals:
