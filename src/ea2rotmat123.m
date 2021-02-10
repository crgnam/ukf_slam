function [A] = ea2rotmat123(roll, pitch, yaw)
% This function calculates the shuster quaternion attitude of the input
% roll, pitch yaw euler angles

    T1 = [1       0             0;
          0   cosd(roll)   sind(roll);
          0  -sind(roll)   cosd(roll)];
          
    T2 = [cosd(pitch)  0  -sind(pitch);
              0      1       0;
          sind(pitch)  0   cosd(pitch)];

    T3 = [ cosd(yaw)   sind(yaw)  0;
          -sind(yaw)   cosd(yaw)  0;
               0           0      1];

    A = T3*T2*T1;
end