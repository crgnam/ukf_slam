%% NEED TO CLEAN UP BEFORE USING AGAIN:
matlabrc; clc; close all;

% Read in the .obj file
obj = readObj('Bennu-Radar.obj');
F = obj.f.v;
V = obj.v*1000;

% Generate points:
min_x = min(V(:,1)); min_y = min(V(:,2)); min_z = min(V(:,3));
max_x = max(V(:,1)); max_y = max(V(:,2)); max_z = max(V(:,3));
num = 50;
[X,Y,Z] = meshgrid(linspace(min_x,max_x,num),linspace(min_y,max_y,num), linspace(min_z,max_z,num));
pts = [X(:),Y(:),Z(:)];

% Isolate points inside the body:
fv.faces = F;
fv.vertices = V;
IN = inpolyhedron(fv, pts);
num_particles = length(pts(IN,1));
disp(num_particles)

% Calculate the gravitational field:
rho = 300;
num_angles = 100;
az = linspace(0,2*pi,num_angles);
el = linspace(-pi/2,pi/2,num_angles);
accel = zeros(length(el),length(az),3);
accel_mag = zeros(length(el),length(az));

% Total mass of bennu:
M = 1.75*78*1e9; %(kg) Mass scaled so surface gravity is equivalent... weird hack, ignore it lmao
G = 6.674*1e-11; %(N m^2 /kg^2) Gravitational Constant
m = M/num_particles; %(kg) mass per particle
gm = M*G;

% Evaluate the field:
for ii = 1:length(el)
    for jj = 1:length(az)
        [x,y,z] = sph2cart(az(jj),el(ii),rho);
        r = [x,y,z];
        
        % Calculate relative position to all test masses:
        r_rel = pts(IN,:) - r;
        [r_rel_u, r_rel_mag] = normr(r_rel);
        
        % Acceleration due to gravity:
        accel_vec = G*m*(sum(r_rel_u./(r_rel_mag.^2)))';
        
        accel(ii,jj,:) = accel_vec;
        accel_mag(ii,jj) = norm(accel_vec);
    end
end

% Show the field map:
surf(az,el,accel_mag,'EdgeColor','none')
axis equal
view([0 90])
title('Bennu Gravity Field From  Volume Filled Finite Point Cloud')

radius = mean([abs(min(V)), max(V)]);
save('accel_map','accel','az','el','rho','gm','radius')

%% Plot:
% patch('Faces',F,'Vertices',V,'FaceColor','none'); hold on
% plot3(pts(IN,1),pts(IN,2),pts(IN,3),'.g','MarkerSize',10)
% plot3(pts(~IN,1),pts(~IN,2),pts(~IN,3),'xr','MarkerSize',10)
% rotate3d on 
% axis equal
% grid on