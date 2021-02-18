%% NEED TO CLEAN UP BEFORE USING AGAIN
matlabrc; clc; close all;

% Load in the generated truth field:
load('accel_map.mat')

num_degrees = 10; %max_degree;
[Cnm_vec, Snm_vec] = coeffs2vec(zeros(num_degrees+1), zeros(num_degrees+1));
num_coeffs = length(Cnm_vec)+length(Snm_vec);
x0 = zeros(num_coeffs, 1);

my_field = GravityField(270, gm, nan, nan);

% Optimize spherical harmonic model:
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
options = optimoptions('fmincon','Display','iter','FiniteDifferenceStepSize', 1e-10,'OptimalityTolerance',1e-16,...
                       'MaxFunctionEvaluations',10000);
x_out = fmincon(@(x) cost(x,my_field,accel,az,el,rho), x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%% Show Results:
[Cnm,Snm] = vec2coeffs(x_out);

% Evaluate the current field:
accel_out = zeros(size(accel));
for ii = 1:length(el)
    for jj = 1:length(az)
        [x,y,z] = sph2cart(az(jj),el(ii),rho);
        r = [x,y,z]';
        accel_out(ii,jj,:) = my_field.acceleration(r, eye(3), Cnm, Snm);
    end
end

% Show the Results:
subplot(1,3,1)
    [~,n2] = normw(accel);
    surf(az,el,n2,'EdgeColor','none')
    colorbar
    caxis([0.9 1.2]*1e-4)
    axis equal
    title('truth field')
    view([0 90])
subplot(1,3,2)
    [~,n] = normw(accel_out);
    surf(az,el,n,'EdgeColor','none')
    colorbar
    caxis([0.9 1.2]*1e-4)
    axis equal
    view([0 90])
    title('Spherical harmonic representation')
subplot(1,3,3)
    surf(az,el,n2-n,'EdgeColor','none')
    colorbar
    caxis([0.01 1]*1e-5)
    axis equal
    view([0 90])
    title('Error')

% Uncomment these lines to save (if it worked well)
% save('Cnm','Cnm')
% save('Snm','Snm')

    
%% Cost Function
% Definition of cost function for optimization
function [eval] = cost(x,my_field,accel_truth,az,el,rho)
    % Unpack spherical harmonic components:
    [Cnm,Snm] = vec2coeffs(x);
    
    % Evaluate the current field:
    accel = zeros(size(accel_truth));
    for ii = 1:length(el)
        for jj = 1:length(az)
            [x,y,z] = sph2cart(az(jj),el(ii),rho);
            r = [x,y,z]';
            accel(ii,jj,:) = my_field.acceleration(r, eye(3), Cnm, Snm);
        end
    end
    
    % Calculate cost:
    [~,n] = normw(accel*1e6 - accel_truth*1e6);
    eval = sum(sum(n));
%     disp(eval)
end
