clear; matlabrc; clc; close all;
addpath(genpath('../src'))
addpath(genpath('../lib'))

A = randn(3,30);

R = ea2rotmat123(5,3,10);
t = [1; 2; 3]*.1;
scale = .76;

B = scale*R*(A - t);

% Use Kabsch algorithm to estimate transformation between the two:
% [R2, t2, S, lrms] = kabsch(A, B);
[regParams,Bfit,ErrorStats] = absor(A,B,'doScale',true,'doTrans',true);
R2 = regParams.R;
t2 = regParams.t;
S2 = regParams.s;
C = (1/S2)*R2'*(B - t2);

%% Show Results:
% figure(1)
% subplot(1,2,1)
%     plot3(A(1,:),A(2,:),A(3,:),'.b','MarkerSize',20); hold on
%     plot3(B(1,:),B(2,:),B(3,:),'xr','MarkerSize',20)
%     axis equal; grid on
%     
% subplot(1,2,2)
%     plot3(A(1,:),A(2,:),A(3,:),'.b','MarkerSize',20); hold on
%     plot3(C(1,:),C(2,:),C(3,:),'xr','MarkerSize',20)
%     axis equal; grid on
    
figure(2)
    plot3(A(1,:),A(2,:),A(3,:),'.b','MarkerSize',20); hold on
%     plot3(B(1,:),B(2,:),B(3,:),'xr','MarkerSize',20)
    plot3(C(1,:),C(2,:),C(3,:),'og','MarkerSize',20)
    axis equal; grid on
