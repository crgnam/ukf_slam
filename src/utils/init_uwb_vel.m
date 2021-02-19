function [v_uwb] = init_uwb_vel(r_uwb,v)
    % Get rotation for -z pointing nadir:
    rotmat = nadir(r_uwb(:,1),v);
    
    % Get the magnitude of the total orbit:
    num_uwb = size(r_uwb,2);
    v_norm = norm(v);
    
    % Generate directions that are random and perpencidular to nadir:
    vel_rand = [randn(2,num_uwb); zeros(1,num_uwb)];
    vel_rand = rotmat'*vel_rand;
    
    % Calculate the new velocity directions:
    v_uwb = (v_norm*rand_ab(1,num_uwb, 0.9, 1.3)).*normc(v + vel_rand);
end

function [ret] = rand_ab(n,m,a,b)
    ret = a + (b-a)*rand(n,m);
end