function val = sideDispx(s,pp,edgN,arcLen,hScale)
%This function is called by RECTMESH.

sz = size(s);
s  = reshape(s,1,prod(sz));

if edgN == 4 
   %Displacement if the difference between the target and the domain.
   val = ppval(pp,1-s/arcLen)-(1-s/arcLen)*hScale;
elseif edgN == 3 
   val = ppval(pp,1-s/arcLen)-hScale;
elseif edgN == 2 
   val = ppval(pp,s/arcLen)-s/arcLen*hScale;
elseif edgN == 1 
   val = ppval(pp,s/arcLen);
end

val = reshape(val,sz);
