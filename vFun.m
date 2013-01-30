function [f,df] = vFun(v,pp,dpp1,dpp2)

f  = -fnval(pp,v.');
df = -[fnval(dpp1,v.') fnval(dpp2,v.')];