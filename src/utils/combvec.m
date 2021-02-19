function [vec_pairs] = combvec(r1,r2)
    r1_row = r1';
    r2_row = r2';
    r1_rep = repelem(r1_row,size(r2_row,1),1);
    r2_rep = repmat(r2_row,size(r1_row,1),1);
    vec_pairs = [r1_rep,r2_rep];
end