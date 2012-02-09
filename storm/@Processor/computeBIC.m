function bic = computeBIC(logL,nParam,nPoints)

bic = -2*logL+nParam*log(nPoints);

end

