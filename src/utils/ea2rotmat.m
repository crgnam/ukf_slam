function [rotmat] = ea2rotmat(rot_x,rot_y,rot_z,sequence)
    % Parse input:
    if nargin == 3
        sequence = '321';
    end
    
    % Generate basic rotaiton matrices:
    T = zeros(3,3,3);
    T(:,:,1) = [1       0             0;
                0   cosd(rot_x)   sind(rot_x);
                0  -sind(rot_x)   cosd(rot_x)];
    T(:,:,2) = [cosd(rot_y)  0  -sind(rot_y);
                      0      1       0;
                sind(rot_y)  0   cosd(rot_y)];
    T(:,:,3) = [ cosd(rot_z)   sind(rot_z)  0;
                -sind(rot_z)   cosd(rot_z)  0;
                     0           0      1];

    % Calculate overall rotation matrix:
    rotmat = T(:,:,str2double(sequence(3)))*...
             T(:,:,str2double(sequence(2)))*...
             T(:,:,str2double(sequence(1)));
end