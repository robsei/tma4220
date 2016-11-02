function p = getVarigrid_recursion(U, idx1, idx2, minvar, cvc, depth, alldepth, maxdepth)

if depth == maxdepth
    warning('Reached maximal iteration depth.');
    return
end

% Get portion of U where to work on.
U_portion = U(idx1, idx2);
nx1 = length(idx1);
nx2 = length(idx2);

% Variance criterion.
v = var(U_portion(:));

% Corner value criterion.
tl = U(idx1(1), idx2(1));
tr = U(idx1(1), idx2(end));
bl = U(idx1(end), idx2(1));
br = U(idx1(end), idx2(end));

d = abs([tl - tr; tr - br; br - bl; bl - tl; tl - br; tr - bl]);

% Make split decicion (include alldepth).
split = (v > minvar) || (depth <= alldepth) || any(d > cvc);

p = [];
if split && (nx1 >= 2) && (nx2 >= 2)
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
    p = vertcat(p, [idx1_mid idx2_mid]);
    
    p = vertcat(p, getVarigrid_recursion(U, idx1_left, idx2_left, minvar, cvc, depth+1, alldepth, maxdepth));
    p = vertcat(p, getVarigrid_recursion(U, idx1_left, idx2_right, minvar, cvc, depth+1, alldepth, maxdepth));
    p = vertcat(p, getVarigrid_recursion(U, idx1_right, idx2_left, minvar, cvc, depth+1, alldepth, maxdepth));
    p = vertcat(p, getVarigrid_recursion(U, idx1_right, idx2_right, minvar, cvc, depth+1, alldepth, maxdepth));
end

