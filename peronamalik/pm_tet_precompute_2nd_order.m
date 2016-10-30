clear all; clc; close all;

syms a0 a1 a2 a3 a4 a5 a6 a7 a8
syms dx1 dx2 i j x1 x2
syms nx1
syms u_iso_bd(eid) g_ani(eid)

assume([a0 a1 a2 a3 a4 a5 a6 a7 a8 dx1 dx2 i j x1 x2 nx1], 'real')

sub2ind = @(sub) 2 * sub(1) + 2 * nx1 * 2 * sub(2) + 1;
rect_idx = [i j; i j+1/2; i j+1; i+1/2 j; i+1/2 j+1/2; i+1/2 j+1; i+1 j; i+1 j+1/2; i+1 j+1];

for k = 1:9
    for p = 1:9
        eqn(p) = a0 + a1*x1 + a2*x1^2 + a3*x2 + a4*x1*x2 + a5*x1^2*x2 + a6*x2^2 + a7*x1*x2^2 + a8*x1^2*x2^2 == (k == p);
    end

    coeffs{k} = solve(eqn, [a0, a1, a2, a3, a4, a5, a6, a7, a8]);
end

bilin = @(coeff) [coeff.a0 + coeff.a1*x1 + coeff.a2*x1^2 + coeff.a3*x2 + coeff.a4*x1*x2 + coeff.a5*x1^2*x2 + coeff.a6*x2^2 + coeff.a7*x1*x2^2 + coeff.a8*x1^2*x2^2];
deriv = @(coeff) [coeff.a1 + 2*coeff.a2*x1 + coeff.a4*x2 + 2*coeff.a5*x1*x2 + coeff.a7*x2^2 + 2*coeff.a8*x1*x2^2; coeff.a3 + coeff.a4*x1 + coeff.a5*x1^2 + 2*coeff.a6*x2 + 2*coeff.a7*x1*x2 + 2*coeff.a8*x1^2*x2];

for k1 = 1:9
    for k2 = 1:9
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = dot(deriv(coeffs{k1}),deriv(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['A_iso(' char(nid1) ',' char(nid2) ') = A_iso(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end

for k1 = 1:9
    for k2 = 1:9
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = dot(bilin(coeffs{k1}),bilin(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['M(' char(nid1) ',' char(nid2) ') = M(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end

disp(' ');

for k1 = 1:9
    eid = sub2ind(rect_idx(1,:));
    nid = sub2ind(rect_idx(k1,:));
    term = u_iso_bd(nid) * deriv(coeffs{k1});
    term = subs(term, x1, (i+1/2)*dx1);
    term = subs(term, x2, (j+1/2)*dx2);
    disp(['d(' char(eid) ',1) = d(' char(eid) ',1) + ' char(term(1)) ';'])
    disp(['d(' char(eid) ',2) = d(' char(eid) ',2) + ' char(term(2)) ';'])
end

disp(' ');

for k1 = 1:9
    for k2 = 1:9
        eid = sub2ind(rect_idx(1,:));
        nid1 = sub2ind(rect_idx(k1,:));
        nid2 = sub2ind(rect_idx(k2,:));
        term = g_ani(eid) * dot(deriv(coeffs{k1}),deriv(coeffs{k2}));
        term = subs(term, x1, (i+1/2)*dx1);
        term = subs(term, x2, (j+1/2)*dx2);
        disp(['A_ani(' char(nid1) ',' char(nid2) ') = A_ani(' char(nid1) ',' char(nid2) ') + ' char(term * dx1 * dx2) ';'])
    end
end