function u = bndDisp(s,pp)

sz = size(s);

u = ppval(pp,reshape(s,1,prod(sz)));
u = reshape(u,sz);
