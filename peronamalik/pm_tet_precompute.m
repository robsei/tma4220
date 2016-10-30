clear all; clc; close all;

syms a b c d
syms dx1 dx2 i j x1 x2
syms nx1
syms u_iso_bd(eid) g_ani(eid)

assume([a b c d dx1 dx2 i j x1 x2 nx1], 'real')

sub2ind = @(sub) sub(1) + nx1 * sub(2) + 1;
rect_idx = [i j; i j+1; i+1 j+1; i+1 j];

for k = 1:4
    for p = 1:4
        eqn(p) = a*rect_idx(p,1)*dx1 + b*rect_idx(p,2)*dx2 + c*rect_idx(p,1)*rect_idx(p,2)*dx1*dx2 + d == (k == p);
    end

    coeffs{k} = solve(eqn, [a, b, c, d]);
end

bilin = @(coeff) [coeff.a * x1 + coeff.b * x2 + coeff.c * x1 * x2 + coeff.d];
deriv = @(coeff) [coeff.a + coeff.c * x2; coeff.b + coeff.c * x1];

for k1 = 1:4
    for k2 = 1:4
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = dot(deriv(coeffs{k1}),deriv(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['A_iso(' char(nid1) ',' char(nid2) ') = A_iso(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end

for k1 = 1:4
    for k2 = 1:4
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = dot(bilin(coeffs{k1}),bilin(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['M(' char(nid1) ',' char(nid2) ') = M(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end

disp(' ');

for k1 = 1:4
    eid = sub2ind(rect_idx(1,:));
    nid = sub2ind(rect_idx(k1,:));
    term = u_iso_bd(nid) * deriv(coeffs{k1});
    term = subs(term, x1, (i+1/2)*dx1);
    term = subs(term, x2, (j+1/2)*dx2);
    disp(['d(' char(eid) ',1) = d(' char(eid) ',1) + ' char(term(1)) ';'])
    disp(['d(' char(eid) ',2) = d(' char(eid) ',2) + ' char(term(2)) ';'])
end

disp(' ');

for k1 = 1:4
    for k2 = 1:4
        eid = sub2ind(rect_idx(1,:));
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = g_ani(eid) * dot(deriv(coeffs{k1}),deriv(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['A_ani(' char(nid1) ',' char(nid2) ') = A_ani(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end