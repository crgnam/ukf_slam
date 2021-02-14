function [X_hat] = initializeMap(camera,image_lmks,rotmat,r_hat,radius_estimate)
    % This function is used to produce initial map estimates by proejcting
    % their image coordinates back onto a sphere:

    % Generate the rays:
    rays    = camera.generateRays(image_lmks,rotmat);
    rays = rays';
    origins = repmat(r_hat',size(rays,1),1);

    % Trace the rays to find ray-sphere intersections:
    lines  = [origins rays];
    sphere = [0 0 0 radius_estimate];
    intersects = intersectLineSphere(lines, sphere);

    % Get the closest ray intersections for each set:
    intersects = getClosestIntersects(intersects,r_hat);
    
    % Store as a new state vector:
    X_hat = [intersects(:,1); intersects(:,2); intersects(:,3)];
end