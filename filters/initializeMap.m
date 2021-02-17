function [X_hat] = initializeMap(camera,image_lmks,rotmat,r_hat,radius_estimate,radius_max)
    % This function is used to produce initial map estimates by proejcting
    % their image coordinates back onto a sphere:

    % Calculate distance:
    distance = norm(r_hat);
    
    % Generate the rays:
    rays = camera.generateRays(image_lmks,rotmat);
    rays = rays';
    origins = zeros(size(rays,1),3);

    % Trace the rays to find ray-sphere intersections:
    lines  = [origins rays];
    sphere_pos = rotmat'*[0;0;-distance];
    sphere = [sphere_pos' radius_estimate];
    intersects = intersectLineSphere(lines, sphere);
    
    % If the sphere intersection failed, do a simply distance projection
    if any(isnan(intersects))
        sphere = [sphere_pos' (radius_max+radius_estimate)/2];
        intersects = intersectLineSphere(lines, sphere);
%         projectedPoints = origins + distance*rays - sphere_pos';
%         X_hat = [projectedPoints(:,1);
%                  projectedPoints(:,2);
%                  projectedPoints(:,3)];
%         disp('SPHERE PROJECTION FAILED.  USING SIMPLE RAY PROJECTION')
        disp('SPHERE PROJECTION FAILED.  INCREASING SIZE OF SPHERE')
    end
    % Get the closest ray intersections for each set:
    intersects = intersects - sphere_pos';
    intersects = getClosestIntersects(intersects,r_hat);

    % Store as a new state vector:
    X_hat = [intersects(:,1); intersects(:,2); intersects(:,3)];
end