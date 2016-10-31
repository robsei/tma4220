clear all; clc; close all;

syms a b c d
syms i j x1 x2
syms nx1
syms u_iso_bd(eid) g_ani(eid)

assume([a b c d i j x1 x2 nx1], 'real')

sub2ind = @(i,j) i + nx1 * (j-1);
rect_idx = [i j; i j+1; i+1 j; i+1 j+1];

% Have four basis functions (four nodes) in linear system.
for k = 1:4
    % Have four contraints in linear system.
    for p = 1:4
        eqn(p) = a*rect_idx(p,1) + b*rect_idx(p,2) + c*rect_idx(p,1)*rect_idx(p,2) + d == (k == p);
    end

    coeffs{k} = solve(eqn, [a, b, c, d]);
end

bilin = @(coeff) coeff.a * x1 + coeff.b * x2 + coeff.c * x1 * x2 + coeff.d;
deriv = @(coeff) [coeff.a + coeff.c * x2; coeff.b + coeff.c * x1];

for n1 = 1:4
    for n2 = 1:4
        nid1 = sub2ind(rect_idx(n1,1),rect_idx(n1,2));
        nid2 = sub2ind(rect_idx(n2,1),rect_idx(n2,2));
        term = dot(deriv(coeffs{n1}),deriv(coeffs{n2}));
        term = subs(term, x1, i+1/2);
        term = subs(term, x2, j+1/2);
        disp(['A_iso(' char(nid1) ',' char(nid2) ') = A_iso(' char(nid1) ',' char(nid2) ') + ' char(term) ';'])
    end
end

for n1 = 1:4
    for n2 = 1:4
        nid1 = sub2ind(rect_idx(n1,1),rect_idx(n1,2));
        nid2 = sub2ind(rect_idx(n2,1),rect_idx(n2,2));
        term = dot(bilin(coeffs{n1}),bilin(coeffs{n2}));
        term = subs(term, x1, i+1/2);
        term = subs(term, x2, j+1/2);
        disp(['M(' char(nid1) ',' char(nid2) ') = M(' char(nid1) ',' char(nid2) ') + ' char(term) ';'])
    end
end

disp(' ');

for n = 1:4
    eid = sub2ind(rect_idx(1,1),rect_idx(1,2));
    nid = sub2ind(rect_idx(n,1),rect_idx(n,2));
    term = u_iso_bd(nid) * deriv(coeffs{n});
    term = subs(term, x1, i+1/2);
    term = subs(term, x2, j+1/2);
    disp(['d(' char(eid) ',1) = d(' char(eid) ',1) + ' char(term(1)) ';'])
    disp(['d(' char(eid) ',2) = d(' char(eid) ',2) + ' char(term(2)) ';'])
end

disp(' ');

for n1 = 1:4
    for n2 = 1:4
        eid = sub2ind(rect_idx(1,1),rect_idx(1,2));
        nid1 = sub2ind(rect_idx(n1,1),rect_idx(n1,2));
        nid2 = sub2ind(rect_idx(n2,1),rect_idx(n2,2));
        term = g_ani(eid) * dot(deriv(coeffs{n1}),deriv(coeffs{n2}));
        term = subs(term, x1, i+1/2);
        term = subs(term, x2, j+1/2);
        disp(['A_ani(' char(nid1) ',' char(nid2) ') = A_ani(' char(nid1) ',' char(nid2) ') + ' char(term) ';'])
    end
end