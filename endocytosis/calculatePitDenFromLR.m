function [currDen] = calculatePitDenFromLR(currKR,distance)
[dlx,dly] = size(currKR);
carea = distance.^2;
careadiff = carea; careadiff(2:length(careadiff)) = diff(carea);
amat = repmat(careadiff',1,dly);
currKRdiff = currKR;
currKRdiff(2:length(distance),:) = diff(currKR,1);
currDen = currKRdiff./amat;
end %of function calculatePitDenFromLR
