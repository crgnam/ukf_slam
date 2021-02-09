function [A] = e2a(roll, pitch, yaw)
% This function calculates the shuster quaternion attitude of the input
% roll, pitch yaw euler angles

%     if roll == 180
%         roll = 179.99;
%     end
%     if pitch == 180
%         pitch = 179.99;
%     end
%     if yaw == 180
%         yaw = 179.99;
%     end

    T1 = [1       0             0;
          0   cosd(roll)   sind(roll);
          0  -sind(roll)   cosd(roll)];
          
    T2 = [cosd(pitch)  0  -sind(pitch);
              0      1       0;
          sind(pitch)  0   cosd(pitch)];

    T3 = [ cosd(yaw)   sind(yaw)  0;
          -sind(yaw)   cosd(yaw)  0;
               0           0      1];

    A = T1*T2*T3;
end