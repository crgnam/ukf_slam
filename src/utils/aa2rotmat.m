function [rotmat] = aa2rotmat(axis, theta)
    rotmat = cos(theta)*eye(3) + (1 - cos(theta))*(axis*axis') - sin(theta)*cpm(axis);
end