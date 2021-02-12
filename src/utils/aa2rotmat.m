function [R] = aa2rotmat(axis, theta)
    % Make sure that the axis is a unit vector
    axis = axis/norm(axis);

    % Construct the rotation matrix
    R = cosd(theta)*eye(3) + (1 - cosd(theta))*(axis*axis') - sind(theta)*cpm(axis);
end