function checkVector=cutIntfWithMask(mask,curve)
max_x=max(max(mask(:,1)),max(curve(:,1)))+1;
max_y=max(max(mask(:,2)),max(curve(:,2)))+1;
mask(max_y,max_x)=0;
indCurve=sub2ind(size(mask),curve(:,2), curve(:,1));
checkVector=mask(indCurve);
%indHits = find(mask(indCurve));

return;


