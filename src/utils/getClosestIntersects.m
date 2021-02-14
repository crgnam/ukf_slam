function [intersects_out] = getClosestIntersects(intersects,cam_origin)
    % Separate out each pair of intersections:
    N = size(intersects,1);
    ints1 = intersects(1:N/2,:);
    ints2 = intersects(N/2+1:end,:);
    
    % Check the lengths of each:
    [~,n1] = normr(ints1-cam_origin');
    [~,n2] = normr(ints2-cam_origin');
    
    intersects_out = zeros(size(ints1));
    inds = n1<n2;
    intersects_out(inds,:)  = ints1(inds,:);
    intersects_out(~inds,:) = ints2(~inds,:);
end