function u = maskApply(U, p);

u = U(sub2ind(size(U), p(:,1), p(:,2)));

end
