function val = sideDispy(s,pp,edgN,arcLen,vScale)
%This function is called by RECTMESH.

sz = size(s);
s  = reshape(s,1,prod(sz));

if edgN == 4 
   val = ppval(pp,1-s/arcLen);
elseif edgN == 3 
   val = ppval(pp,1-s/arcLen)-(1-s/arcLen)*vScale;
elseif edgN == 2 
   val = ppval(pp,s/arcLen)-vScale;
elseif edgN == 1 
   val = ppval(pp,s/arcLen)-s/arcLen*vScale;
end

val = reshape(val,sz);
