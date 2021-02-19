function [vec_pairs] = combvec(r)
    r_row = r';
    num_rows = size(r_row,1);
    
    % Generate all pairs of vectors:
    [A,B] = meshgrid(1:num_rows,1:num_rows);
    c = cat(2,A',B');
    d = reshape(c,[],2);
    
    % Remove pairs that are of the same vector:
    d(d(:,1) == d(:,2),:) = [];
    
    % Create vector pairs:
    vec_pairs = [r_row(d(:,1),:),r_row(d(:,2),:)];
end