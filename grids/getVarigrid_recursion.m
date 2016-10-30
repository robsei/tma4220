function p_idx = getVarigrid_recursion(U, idx1, idx2, minvar, depth, maxdepth)

u = U(idx1, idx2);
p_idx = [];

nx1 = length(idx1);
nx2 = length(idx2);

if depth == maxdepth
    warning('Reached maximal iteration depth.');
    return
end

v = var(u(:));

if (v > minvar) && (nx1 >= 2) && (nx2 >= 2)
    % Split indices into halfs.   
    idx1_split = floor(nx1/2);
    idx1_left = idx1(1:idx1_split);
    idx1_right = idx1(idx1_split+1:end);
    idx1_mid = idx1_right(1);

    idx2_split = floor(nx2/2);
    idx2_left = idx2(1:idx2_split);
    idx2_right = idx2(idx2_split+1:end);
    idx2_mid = idx2_right(1);
    
    % Add midpoint of variance-rich area to list of points.
    p_idx = vertcat(p_idx, [idx1_mid idx2_mid]);
    
    p_idx = vertcat(p_idx, getVarigrid_recursion(U, idx1_left, idx2_left, minvar, depth+1, maxdepth));
    p_idx = vertcat(p_idx, getVarigrid_recursion(U, idx1_left, idx2_right, minvar, depth+1, maxdepth));
    p_idx = vertcat(p_idx, getVarigrid_recursion(U, idx1_right, idx2_left, minvar, depth+1, maxdepth));
    p_idx = vertcat(p_idx, getVarigrid_recursion(U, idx1_right, idx2_right, minvar, depth+1, maxdepth));
end

