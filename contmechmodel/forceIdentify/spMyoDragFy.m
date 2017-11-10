function fy = spMyoDragFy(t,x,y,sp)

sz = size(x);

fy = reshape(fnval(sp,[reshape(x,1,prod(sz));reshape(y,1,prod(sz))]));
