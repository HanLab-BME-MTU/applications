function fx = spMyoDragFx(t,x,y,sp)

sz = size(x);

fx = reshape(fnval(sp,[reshape(x,1,prod(sz));reshape(y,1,prod(sz))]),sz);

